import logging
import math
import os
import sys

from configobj import ConfigObj
from openmm.unit import kelvin, kilocalories_per_mole

from gibbs_sampling import pairwise_independence_sampling
from local_openmm_transport_sync import LocalOpenMMTransport
from ommreplica import OMMReplicaATM
from ommsystem import OMMSystemAmberRBFE
from ommworker_sync import OMMWorkerATM


class sync_re:
    """
    Class to set up and run asynchronous file-based RE calculations
    """
    logging.config.fileConfig(os.path.join(os.path.dirname(__file__), "utils/logging.conf"))

    def __init__(self, config_file, options):
        self._setLogger()

        self.command_file = config_file

        self.jobname = os.path.splitext(os.path.basename(config_file))[0]
        self.config = ConfigObj(self.command_file)

        self._checkInput()
        self._printStatus()

    def _setLogger(self):
        self.logger = logging.getLogger("async_re")

    def __getattribute__(self, name):
        if name == 'replicas_waiting':
            # Return a list of replica indices of replicas in a wait state.
            return [k for k in range(self.nreplicas)
                    if self.status[k]['running_status'] == 'W']
        elif name == 'states_waiting':
            # Return a list of state ids of replicas in a wait state.
            return [self.status[k]['stateid_current']
                    for k in self.replicas_waiting]
        elif name == 'replicas_waiting_to_exchange':
            # Return a list of replica indices of replicas in a wait state that
            # have ALSO completed at least one cycle.
            return [k for k in range(self.nreplicas)
                    if (self.status[k]['running_status'] == 'W' and
                        self.status[k]['cycle_current'] > 1)]
        elif name == 'states_waiting_to_exchange':
            # Return a list of state ids of replicas in a wait state that have
            # ALSO completed at least one cycle.
            return [self.status[k]['stateid_current']
                    for k in self.replicas_waiting_to_exchange]
        elif name == 'waiting':
            return len(self.replicas_waiting)
        elif name == 'replicas_running':
            # Return a list of replica indices of replicas in a running state.
            return [k for k in range(self.nreplicas)
                    if self.status[k]['running_status'] == 'R']
        elif name == 'running':
            return len(self.replicas_running)
        else:
            return object.__getattribute__(self,name)

    def _printStatus(self):
        """Print a report of the input parameters."""
        self.logger.info("ASyncRE-OpenMM, Version")
        self.logger.info("command_file = %s", self.command_file)
        self.logger.info("jobname = %s", self.jobname)
        self.logger.info("Keywords:")
        for k,v in self.config.iteritems():
            self.logger.info("%s: %s", k, v)

    def _checkInput(self):
        """
        Check that required parameters are specified. Parse these and other
        optional settings.
        """
        # Required Options
        #
        # basename for the job
        self.basename = self.config.get('BASENAME')
        assert self.basename, 'BASENAME needs to be specified'

        #job transport mechanism
        self.transport_mechanism = self.config.get('JOB_TRANSPORT')
        self.transport = None

        if self.transport_mechanism == "LOCAL_OPENMM":

            nodefile = self.config.get('NODEFILE')
            assert nodefile, "NODEFILE needs to be specified"
            """
            check the information in the nodefile. there should be six columns in the  
            nodefile. They are 'node name', 'slot number', 'number of threads', 
            'platform','username', and 'name of the temperary folder'
            """
            node_info= []
            with open(nodefile, 'r') as f:
                line=f.readline()
                nodeid = 0
                while line:
                    lineID=line.split(",")
                    node_info.append({})
                    node_info[nodeid]["node_name"] = str(lineID[0].strip())
                    node_info[nodeid]["slot_number"] = str(lineID[1].strip())
                    node_info[nodeid]["threads_number"] = str(lineID[2].strip())
                    node_info[nodeid]["arch"] = str(lineID[3].strip())
                    node_info[nodeid]["user_name"] = str(lineID[4].strip())
                    node_info[nodeid]["tmp_folder"] = str(lineID[5].strip())
                    #tmp_folder has to be pre-assigned
                    assert node_info[nodeid]["tmp_folder"] != "", 'tmp_folder in nodefile needs to be specified'
                    nodeid += 1
                    line = f.readline()

            #set the nodes information
            self.num_nodes = len(node_info)
            self.compute_nodes = node_info
            #Can print out here to check the node information
            self.logger.info("compute nodes: %s", ', '.join([n['node_name'] for n in node_info]))

        # number of replicas (may be determined by other means)
        self.nreplicas = None

        # verbose printing
        if self.config.get('VERBOSE').lower() == 'yes':
            self.verbose = True
            if self.logger:
                self.logger.setLevel(logging.DEBUG)
        else:
            self.verbose = False

        # self.jobname = self.config.get('BASENAME')

    def setupJob(self):
        self.transport = LocalOpenMMTransport(self.basename, self.openmm_workers, self.openmm_replicas)
        # create status table
        self.status = [{'stateid_current': k, 'running_status': 'W',
                        'cycle_current': 1} for k in range(self.nreplicas)]
        for replica in self.openmm_replicas:
            self.status[replica._id]['cycle_current'] = replica.get_cycle()
            self.status[replica._id]['stateid_current'] = replica.get_stateid()
            self.logger.info("Replica %d Cycle %d Stateid %d" % (replica._id, self.status[replica._id]['cycle_current'], self.status[replica._id]['stateid_current']))

        self.updateStatus()
        self.print_status()

    def scheduleJobs(self):

        max_samples = None
        assert self.config.get('MAX_SAMPLES'), "MAX_SAMPLES has to be specified"
        max_samples = int(self.config.get('MAX_SAMPLES'))

        # TODO: restart properly
        num_sims = max_samples * len(self.openmm_replicas)
        for i_sim in range(num_sims):
            self.logger.info(f"Simulation: {i_sim+1}/{num_sims}")

            self.updateStatus()
            self.print_status()
            self.launchJobs()
            self.updateStatus()
            self.print_status()

            self.updateStatus()
            self.print_status()
            self.doExchanges()
            self.print_status()

            self.checkpointJob()

    def updateStatus(self):
        """Scan the replicas and update their states."""
        # self.transport.poll()
        for k in range(self.nreplicas):
            # self._updateStatus_replica(k)
            self.update_state_of_replica(k)

    def _cycle_of_replica(self,repl):
        return self.status[repl]['cycle_current']

    def launchJobs(self):
        wait = sorted(self.replicas_waiting, key=self._cycle_of_replica)
        k = wait[0]
        self.logger.info('Launching replica %d cycle %d', k, self.status[k]['cycle_current'])
        self._launchReplica(k,self.status[k]['cycle_current'])
        self.status[k]['cycle_current'] += 1

    def doExchanges(self):
        self.logger.info("Replica exchange")

        replicas_to_exchange = self.replicas_waiting_to_exchange
        states_to_exchange = self.states_waiting_to_exchange

        self.logger.debug(f"Replicas to exchange: {replicas_to_exchange}")
        self.logger.debug(f"States to exchange: {states_to_exchange}")

        if len(replicas_to_exchange) < 2:
            return 0

        # Matrix of replica energies in each state.
        # The computeSwapMatrix() function is defined by application classes
        swap_matrix = self._computeSwapMatrix(replicas_to_exchange, states_to_exchange)
        self.logger.debug(swap_matrix)

        for repl_i in replicas_to_exchange:
            sid_i = self.status[repl_i]['stateid_current']
            curr_states = [self.status[repl_j]['stateid_current']
                           for repl_j in replicas_to_exchange]
            repl_j = pairwise_independence_sampling(repl_i,sid_i,
                                                    replicas_to_exchange,
                                                    curr_states,
                                                    swap_matrix)
            if repl_j != repl_i:
                sid_i = self.status[repl_i]['stateid_current']
                sid_j = self.status[repl_j]['stateid_current']
                self.status[repl_i]['stateid_current'] = sid_j
                self.status[repl_j]['stateid_current'] = sid_i
                self.logger.info("Replica %d new state %d" % (repl_i, sid_j))
                self.logger.info("Replica %d new state %d" % (repl_j, sid_i))


class openmm_job(sync_re):
    def __init__(self, command_file, options):
        sync_re.__init__(self, command_file, options)
        self.openmm_replicas = None
        self.stateparams = None
        self.openmm_workers = None
        self.kb = 0.0019872041*kilocalories_per_mole/kelvin

    def _setLogger(self):
        self.logger = logging.getLogger("async_re.openmm_sync_re")

    def checkpointJob(self):
        for replica in self.openmm_replicas:
            replica.save_checkpoint()

    def _launchReplica(self,replica,cycle):

        nsteps = int(self.config.get('PRODUCTION_STEPS'))
        nprnt = int(self.config.get('PRNT_FREQUENCY'))
        ntrj = int(self.config.get('TRJ_FREQUENCY'))
        assert nprnt % nsteps == 0, "PRNT_FREQUENCY must be an integer multiple of PRODUCTION_STEPS."
        assert ntrj % nsteps == 0, "TRJ_FREQUENCY must be an integer multiple of PRODUCTION_STEPS."

        job_info = {"replica": replica, "cycle": cycle, "nsteps": nsteps, "nprnt": nprnt, "ntrj": ntrj}

        return self.transport.launchJob(replica, job_info)

    def update_state_of_replica(self, repl):
        replica = self.openmm_replicas[repl]
        #retrieve previous state if set
        (old_stateid, old_par) =  replica.get_state()
        if old_stateid != None:
            old_temperature = old_par['temperature']
        #sets new state
        stateid = self.status[repl]['stateid_current']
        par = self.stateparams[stateid]
        replica.set_state(stateid, par)

        #rescale velocities (relevant only if state has changed)
        if old_stateid != None:
            if stateid != old_stateid:
                temperature = par['temperature']
                scale = math.sqrt(temperature/old_temperature)
                for i in range(0,len(replica.velocities)):
                    replica.velocities[i] = scale*replica.velocities[i]

    def _getPar(self, repl):
        _, par = self.openmm_replicas[repl].get_state()
        return par

    def _computeSwapMatrix(self, repls, states):
        """
        Compute matrix of dimension-less energies: each column is a replica
        and each row is a state so U[i][j] is the energy of replica j in state
        i.

        Note that the matrix is sized to include all of the replicas and states
        but the energies of replicas not in waiting state, or those of waiting
        replicas for states not belonging to waiting replicas list are
        undefined.
        """
        # U will be sparse matrix, but is convenient bc the indices of the
        # rows and columns will always be the same.
        U = [[ 0. for j in range(self.nreplicas)]
             for i in range(self.nreplicas)]

        n = len(repls)

        #collect replica parameters and potentials
        par = []
        pot = []
        for k in repls:
            v = self._getPot(k)
            l = self._getPar(k)
            par.append(l)
            pot.append(v)

        for i in range(n):
            repl_i = repls[i]
            for j in range(n):
                sid_j = states[j]
                #energy of replica i in state j
                U[sid_j][repl_i] = self._reduced_energy(par[j],pot[i])
        return U


class openmm_job_ATM(openmm_job):
    def _buildStates(self):
        self.stateparams = []
        for (lambd,direction,intermediate,lambda1,lambda2,alpha,u0,w0) in zip(self.lambdas,self.directions,self.intermediates,self.lambda1s,self.lambda2s,self.alphas,self.u0s,self.w0coeffs):
            for tempt in self.temperatures:
                par = {}
                par['lambda'] = float(lambd)
                par['atmdirection'] = float(direction)
                par['atmintermediate'] = float(intermediate)
                par['lambda1'] = float(lambda1)
                par['lambda2'] = float(lambda2)
                par['alpha'] = float(alpha)/kilocalories_per_mole
                par['u0'] = float(u0)*kilocalories_per_mole
                par['w0'] = float(w0)*kilocalories_per_mole
                par['temperature'] = float(tempt)*kelvin
                self.stateparams.append(par)
        return len(self.stateparams)

    def _checkInput(self):
        sync_re._checkInput(self)

        assert self.config.get('LAMBDAS'), "LAMBDAS needs to be specified"
        self.lambdas = self.config.get('LAMBDAS').split(',')
        #list of temperatures
        assert self.config.get('TEMPERATURES'), "TEMPERATURES needs to be specified"
        self.temperatures = self.config.get('TEMPERATURES').split(',')

        #flag to identify the intermediate states, typically the one at lambda=1/2
        self.intermediates = self.config.get('INTERMEDIATE').split(',')

        #direction of transformation at each lambda
        #ABFE 1 from RA to R+A, -1 from R+A to A
        #RBFE 1 from RA+B to RB+A, -1 from RB+A to RA+B
        self.directions = self.config.get('DIRECTION').split(',')

        #parameters of the softplus alchemical potential
        #lambda1 = lambda2 gives the linear potential
        self.lambda1s = self.config.get('LAMBDA1').split(',')
        self.lambda2s = self.config.get('LAMBDA2').split(',')
        self.alphas = self.config.get('ALPHA').split(',')
        self.u0s = self.config.get('U0').split(',')
        self.w0coeffs = self.config.get('W0COEFF').split(',')

        #build parameters for the lambda/temperatures combined states
        self.nreplicas = self._buildStates()

    def print_status(self):
        """
        Writes to BASENAME_stat.txt a text version of the status of the RE job
        """
        logfile = "%s_stat.txt" % self.basename
        ofile = open(logfile,"w")
        log = "Replica  State  Lambda Lambda1 Lambda2 Alpha U0 W0coeff Temperature Status  Cycle \n"
        for k in range(self.nreplicas):
            stateid = self.status[k]['stateid_current']
            log += "%6d   %5d  %6.3f %6.3f %6.3f %6.3f %6.2f %6.2f %6.2f %5s  %5d\n" % (k, stateid, self.stateparams[stateid]['lambda'], self.stateparams[stateid]['lambda1'], self.stateparams[stateid]['lambda2'], self.stateparams[stateid]['alpha']*kilocalories_per_mole, self.stateparams[stateid]['u0']/kilocalories_per_mole, self.stateparams[stateid]['w0']/kilocalories_per_mole, self.stateparams[stateid]['temperature']/kelvin, self.status[k]['running_status'], self.status[k]['cycle_current'])
        log += "Running = %d\n" % self.running
        log += "Waiting = %d\n" % self.waiting

        ofile.write(log)
        ofile.close()

    #evaluates the softplus function
    def _softplus(self, lambda1, lambda2, alpha, u0, w0, uf):
        ee = 1.0 + math.exp(-alpha*(uf-u0))
        softplusf = lambda2 * uf + w0
        if alpha._value > 0.:
            softplusf += ((lambda2 - lambda1)/alpha) * math.log(ee)
        return softplusf

    #customized getPot to return the unperturbed potential energy
    #of the replica U0 = U - W_lambda(u)
    def _getPot(self, repl):
        replica = self.openmm_replicas[repl]
        pot = replica.get_energy()
        epot = pot['potential_energy']
        pertpot = pot['perturbation_energy']
        (stateid, par) = replica.get_state()
        # direction = par['atmdirection']
        lambda1 = par['lambda1']
        lambda2 = par['lambda2']
        alpha = par['alpha']
        u0 = par['u0']
        w0 = par['w0']
        ebias = self._softplus(lambda1, lambda2, alpha, u0, w0, pertpot)
        pot['unbiased_potential_energy'] = epot - ebias
        pot['direction'] = par['atmdirection']
        pot['intermediate'] = par['atmintermediate']
        return pot

    def _reduced_energy(self, par, pot):
        temperature = par['temperature']
        beta = 1./(self.kb*temperature)
        direction = par['atmdirection']
        lambda1 = par['lambda1']
        lambda2 = par['lambda2']
        alpha = par['alpha']
        u0 = par['u0']
        w0 = par['w0']
        state_direction = par['atmdirection']
        state_intermediate = par['atmintermediate']
        epot0 = pot['unbiased_potential_energy']
        pertpot = pot['perturbation_energy']
        replica_direction = pot['direction']
        replica_intermediate = pot['intermediate']
        if (replica_direction == state_direction) or (state_intermediate > 0 and replica_intermediate > 0):
            ebias = self._softplus(lambda1, lambda2, alpha, u0, w0, pertpot)
            return beta*(epot0 + ebias)
        else:
            #prevent exchange
            large_energy = 1.e12
            return large_energy


class openmm_job_AmberRBFE(openmm_job_ATM):
    def __init__(self, command_file, options):
        super().__init__(command_file, options)

        prmtopfile = self.basename + ".prmtop"
        crdfile = self.basename + ".inpcrd"

        if self.stateparams is None:
            self._buildStates()

        #builds service worker for replicas use
        service_ommsys = OMMSystemAmberRBFE(self.basename, self.config, prmtopfile, crdfile, self.logger)
        self.service_worker = OMMWorkerATM(self.basename, service_ommsys, self.config, compute = False, logger = self.logger)
        #creates openmm replica objects
        self.openmm_replicas = []
        for i in range(self.nreplicas):
            replica = OMMReplicaATM(i, self.basename, self.service_worker, self.logger)
            if replica.stateid == None:
                replica.set_state(i, self.stateparams[i])#initial setting
            self.openmm_replicas.append(replica)

        # creates openmm context objects
        self.openmm_workers = []
        for node in self.compute_nodes:
            ommsys = OMMSystemAmberRBFE(self.basename, self.config, prmtopfile, crdfile, self.logger) 
            self.openmm_workers.append(OMMWorkerATM(self.basename, ommsys, self.config, node_info = node, compute = True, logger = self.logger))
