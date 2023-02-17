import logging
import math
import multiprocessing as mp
import os
import signal
import sys
import time

from configobj import ConfigObj
from openmm.unit import kelvin, kilocalories_per_mole

from gibbs_sampling import pairwise_independence_sampling
from local_openmm_transport_sync import LocalOpenMMTransport
from ommreplica import OMMReplicaATM
from ommsystem import OMMSystemAmberRBFE
from ommworker_sync import OMMWorkerATM

__version__ = '3.3.0'


class sync_re:
    """
    Class to set up and run asynchronous file-based RE calculations
    """
    logging.config.fileConfig(os.path.join(os.path.dirname(__file__), "utils/logging.conf"))

    def __init__(self, command_file, options):
        try:
            import setproctitle
            setproctitle.setproctitle("AToM %s" % command_file)
        except:
            pass

        self._setLogger()

        self.command_file = command_file
        if not os.path.exists(self.command_file):
           self._exit('No such file: %s'%self.command_file)

        self.jobname = os.path.splitext(os.path.basename(command_file))[0]
        self.keywords = ConfigObj(self.command_file)

        self._checkInput()
        self._printStatus()

        #set to False to run without exchanges
        self.exchange = True

        #catch ctrl-C and SIGTERM to terminate threads gracefully
        signal.signal(signal.SIGINT, self._signal_handler)
        signal.signal(signal.SIGTERM, self._signal_handler)

    def _exit(self, message):
        """Print and flush a message to stdout and then exit."""
        self.checkpointJob()
        self._cleanup()
        self.logger.info(message)
        sys.stdout.flush()
        self.logger.info('Waiting for children to complete ...')
        while True:
            time.sleep(1)
            if not mp.active_children():
                break
        time.sleep(20)
        sys.exit(1)

    def getVersion(self):
        return __version__

    def _cleanup(self):
        try:
            for worker in self.openmm_workers:
                worker.finish()
        except:
            pass

    def _signal_handler(self, sig, frame):
        msg = "Termination signal %s detected ... cleaning up." % str(sig)
        self._exit(msg)

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
        self.logger.info("ASyncRE-OpenMM, Version %s" % self.getVersion())
        self.logger.info("command_file = %s", self.command_file)
        self.logger.info("jobname = %s", self.jobname)
        self.logger.info("Keywords:")
        for k,v in self.keywords.iteritems():
            self.logger.info("%s: %s", k, v)

    def _checkInput(self):
        """
        Check that required parameters are specified. Parse these and other
        optional settings.
        """
        # Required Options
        #
        # basename for the job
        self.basename = self.keywords.get('BASENAME')
        if self.basename is None:
            self._exit('BASENAME needs to be specified')

        #job transport mechanism
        self.transport_mechanism = self.keywords.get('JOB_TRANSPORT')
        if self.transport_mechanism is None:
            self._exit('JOB_TRANSPORT needs to be specified')
        #only LOCAL_OPENMM is supported for now
        if self.transport_mechanism != "LOCAL_OPENMM" :
            self._exit("unknown JOB_TRANSPORT %s" % self.transport_mechanism)
        # reset job transport
        self.transport = None

        if self.transport_mechanism == "LOCAL_OPENMM":
            if self.keywords.get('NODEFILE') is None:
                self._exit("NODEFILE needs to be specified")
            nodefile = self.keywords.get('NODEFILE')
            """
            check the information in the nodefile. there should be six columns in the  
            nodefile. They are 'node name', 'slot number', 'number of threads', 
            'platform','username', and 'name of the temperary folder'
            """
            node_info= []
            try:
                f = open(nodefile, 'r')
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
                    if node_info[nodeid]["tmp_folder"] == "":
                        self._exit('tmp_folder in nodefile needs to be specified')
                    nodeid += 1
                    line = f.readline()
                f.close()
            except:
                self._exit("Unable to process nodefile %s" % nodefile)
                # reset job transport
                self.transport = None
            #set the nodes information
            self.num_nodes = len(node_info)
            self.compute_nodes=node_info
            #Can print out here to check the node information
            self.logger.info("compute nodes: %s", ', '.join([n['node_name'] for n in node_info]))

        # execution time in minutes
        self.walltime = float(self.keywords.get('WALL_TIME'))
        if self.walltime is None:
            self._exit('WALL_TIME (in minutes) needs to be specified')

        # Optional variables
        #
        env = self.keywords.get('ENGINE_ENVIRONMENT')
        if env is not None and env != '':
            self.engine_environment = env.split(',')
        else:
            self.engine_environment = []

        # number of replicas (may be determined by other means)
        self.nreplicas = None

        # verbose printing
        if self.keywords.get('VERBOSE').lower() == 'yes':
            self.verbose = True
            if self.logger:
                self.logger.setLevel(logging.DEBUG)
        else:
            self.verbose = False

        self.implicitsolvent =  self.keywords.get('IMPLICITSOLVENT')
        self.totalsteps = self.keywords.get('PRODUCTION_STEPS')
        self.jobname = self.keywords.get('BASENAME')
        self.stepgap = self.keywords.get('PRNT_FREQUENCY')

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
        # Gets the wall clock time for a replica to complete a cycle
        # If unspecified it is estimated as 1% of job wall clock time
        replica_run_time = self.keywords.get('REPLICA_RUN_TIME')
        if self.keywords.get('REPLICA_RUN_TIME') is None:
            replica_run_time = int(round(self.walltime/100.))
        else:
            replica_run_time = int(self.keywords.get('REPLICA_RUN_TIME'))
        # double it to give time for current running processes
        # and newly submitted processes to complete
        replica_run_time *= 2

        # Time in between cycles in seconds
        # If unspecified it is set as 30 secs
        if self.keywords.get('CYCLE_TIME') is None:
            cycle_time = 30.0
        else:
            cycle_time = float(self.keywords.get('CYCLE_TIME'))

        if self.keywords.get('MIN_TIME') is None:
            min_time = 1
        else:
            min_time = float(self.keywords.get('MIN_TIME'))

        if self.keywords.get('CHECKPOINT_TIME') is None:
            checkpoint_time = cycle_time
        else:
            checkpoint_time = float(self.keywords.get('CHECKPOINT_TIME'))

        sample_steps = int(self.keywords.get('PRNT_FREQUENCY'))
        cycle_steps = int(self.keywords.get('PRODUCTION_STEPS'))
        cycle_to_sample = sample_steps/cycle_steps
        enough_samples = False
        if self.keywords.get('MAX_SAMPLES')  is not None:
            max_samples = int(self.keywords.get('MAX_SAMPLES'))
            enough_samples = all( [ (replica.get_cycle()-1)/cycle_to_sample >=  max_samples
                                    for replica in self.openmm_replicas ]  )
            if enough_samples:
                self.logger.info("All replicas collected the requested number of samples (%d)" % max_samples)

        start_time = time.time()
        end_time = start_time + 60*(self.walltime - replica_run_time)
        last_checkpoint_time = start_time

        self.logger.debug("Entering event loop: 1")
        while ( time.time() < end_time and
                self.transport.numNodesAlive() > 0 and
                not enough_samples ) :

            self.logger.debug("=== Event loop: 1 ===")

            current_time = time.time()

            self.updateStatus()
            self.print_status()
            self.launchJobs()
            self.updateStatus()
            self.print_status()

            self.transport.ProcessJobQueue(min_time,cycle_time)

            self.updateStatus()
            self.print_status()
            if self.exchange:
                self.doExchanges()
            self._write_status()
            self.print_status()

            if current_time - last_checkpoint_time > checkpoint_time:
                self.logger.info("Checkpointing ...")
                self.checkpointJob()
                last_checkpoint_time = current_time
                self.logger.info("done.")

            #terminates if enough samples have been collected
            if self.keywords.get('MAX_SAMPLES')  is not None:
                max_samples = int(self.keywords.get('MAX_SAMPLES'))
                enough_samples = all( [ (replica.get_cycle()-1)/cycle_to_sample >=  max_samples
                                        for replica in self.openmm_replicas ]  )
                if enough_samples:
                    self.logger.info("All replicas collected the requested number of samples (%d)" % max_samples)

        if self.transport.numNodesAlive() <= 0 :
            self.logger.info("No compute devices are alive. Quitting.")
        else:
            if time.time() >= end_time :
                self.logger.info("Requested simulation time completed (%d mins)" % self.walltime)
            self.logger.info("Normal termination.")

        self.transport.DrainJobQueue()
        self.updateStatus()
        self.print_status()
        self.waitJob()
        self.checkpointJob()
        self.cleanJob()

    def waitJob(self):
        # wait until all jobs are complete
        completed = False
        while not completed:
            self.updateStatus()
            completed = True
            for k in range(self.nreplicas):
                if self.status[k]['running_status'] == "R":
                    completed = False
            time.sleep(1)

    def cleanJob(self):
        self._cleanup()
        return

    def checkpointJob(self):
        #defined in subclasses
        pass

    def _write_status(self):
        pass

    def _read_status(self):
        pass

    def print_status(self):
        """
        Writes to BASENAME_stat.txt a text version of the status of the RE job.
        It's fun to follow the progress in real time by doing:
        watch cat BASENAME_stat.txt
        """
        log = 'Replica  State  Status  Cycle \n'
        for k in range(self.nreplicas):
            log += ('%6d   %5d  %5s  %5d \n'%
                    (k,self.status[k]['stateid_current'],
                     self.status[k]['running_status'],
                     self.status[k]['cycle_current']))
        log += 'Running = %d\n'%self.running
        log += 'Waiting = %d\n'%self.waiting

        logfile = '%s_stat.txt'%self.basename
        ofile = _open(logfile,'w')
        ofile.write(log)
        ofile.close()

    def _buildInpFile(self, repl):
        pass

    def updateStatus(self):
        """Scan the replicas and update their states."""
        self.transport.poll()
        for k in range(self.nreplicas):
            self._updateStatus_replica(k)
        self._write_status()

    def _updateStatus_replica(self, replica):
        this_cycle = self.status[replica]['cycle_current']
        if self.status[replica]['running_status'] == 'R':
            if self.transport.isDone(replica,this_cycle):
                self.status[replica]['running_status'] = 'S'
                #MD engine modules implement ways to check for completion.
                #by testing existence of output file, etc.
                if self._hasCompleted(replica,this_cycle):
                    self.status[replica]['cycle_current'] += 1
                else:
                    self.logger.warning('_updateStatus_replica(): restarting replica %s (cycle %s)',
                                        replica, this_cycle)
                self.status[replica]['running_status'] = 'W'
        self.update_state_of_replica(replica)

    def _njobs_to_run(self):
        # size of subjob buffer as a percentage of job slots
        subjobs_buffer_size = self.keywords.get('SUBJOBS_BUFFER_SIZE')
        if subjobs_buffer_size is None:
            subjobs_buffer_size = 0.5
        else:
            subjobs_buffer_size = float(subjobs_buffer_size)

        # launch new replicas if the number of submitted/running subjobs is
        # less than the number of available slots
        # (total_cores/subjob_cores) + 50%
        available_slots = self.num_nodes
        max_njobs_submittable = int((1.+subjobs_buffer_size)*available_slots)
        nlaunch = self.waiting - max(2,self.nreplicas - max_njobs_submittable)
        nlaunch = max(0,nlaunch)
        if self.verbose:
            self.logger.debug('available_slots: %d', available_slots)
            self.logger.debug('max job queue size: %d', max_njobs_submittable)
            self.logger.debug('running/submitted subjobs: %d', self.running)
            self.logger.debug('waiting replicas: %d', self.waiting)
            self.logger.debug('replicas to launch: %d', nlaunch)
        return nlaunch

    def _cycle_of_replica(self,repl):
        return self.status[repl]['cycle_current']

    def launchJobs(self):
        """
        Scans the replicas in wait state and randomly launches them
        """
        jobs_to_launch = self._njobs_to_run()
        if jobs_to_launch > 0:
            #prioritize replicas that are most behind
            wait = sorted(self.replicas_waiting, key=self._cycle_of_replica)
            #  random.shuffle(wait)
            n = min(jobs_to_launch,len(wait))
            for k in wait[0:n]:
                self.logger.info('Launching replica %d cycle %d', k, self.status[k]['cycle_current'])
                # the _launchReplica function is implemented by
                # MD engine modules
                status = self._launchReplica(k,self.status[k]['cycle_current'])
                if status != None:
                    self.status[k]['running_status'] = 'R'

    def doExchanges(self):
        """Perform exchanges among waiting replicas using Gibbs sampling."""

        #check the exchange
        if self.verbose:
            self.logger.info("starting the replica exchange")

        replicas_to_exchange = self.replicas_waiting_to_exchange
        states_to_exchange = self.states_waiting_to_exchange
        nreplicas_to_exchange = len(replicas_to_exchange)
        if nreplicas_to_exchange < 2:
            return 0

        if self.verbose:
            self.logger.debug('Initiating exchanges amongst %d replicas:', nreplicas_to_exchange)

        exchange_start_time = time.time()

        # Matrix of replica energies in each state.
        # The computeSwapMatrix() function is defined by application classes
        matrix_start_time = time.time()
        swap_matrix = self._computeSwapMatrix(replicas_to_exchange,
                                              states_to_exchange)
        if self.verbose:
            self.logger.info(swap_matrix)
        matrix_time = time.time() - matrix_start_time

        sampling_start_time = time.time()

        for repl_i in replicas_to_exchange:
            #repl_i = choice(replicas_to_exchange)
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

        # Uncomment to debug Gibbs sampling:
        # Actual and observed populations of state permutations should match.
        #
        #     self._debug_collect_state_populations(replicas_to_exchange)
        # self._debug_validate_state_populations(replicas_to_exchange,
        #                                        states_to_exchange,U)
        sampling_time = time.time() - sampling_start_time

    # see children classes for specific implementations
    def update_state_of_replica(self, repl):
        pass

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
        #disable ctrl-c
        s = signal.signal(signal.SIGINT, signal.SIG_IGN)
        # update replica objects of waiting replicas
        self.update_replica_states()
        for replica in self.openmm_replicas:
            replica.save_checkpoint()
        signal.signal(signal.SIGINT, s)

    def _launchReplica(self,replica,cycle):
        nsteps = int(self.keywords.get('PRODUCTION_STEPS'))
        nprnt = int(self.keywords.get('PRNT_FREQUENCY'))
        ntrj = int(self.keywords.get('TRJ_FREQUENCY'))
        if nprnt % nsteps != 0:
            self._exit("PRNT_FREQUENCY must be an integer multiple of PRODUCTION_STEPS.")
        if ntrj % nsteps != 0:
            self._exit("TRJ_FREQUENCY must be an integer multiple of PRODUCTION_STEPS.")

        job_info = {
            "replica": replica,
            "cycle": cycle,
            "nsteps": nsteps,
            "nprnt": nprnt,
            "ntrj": ntrj
        }

        status = self.transport.launchJob(replica, job_info)
        return status

    #sync replicas with the current state assignments
    def update_replica_states(self):
        for repl in range(self.nreplicas):
            self.update_state_of_replica(repl)

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

        #additional operations if any
        self._update_state_of_replica_addcustom(replica)

    def _update_state_of_replica_addcustom(self, replica):
        pass

    def _hasCompleted(self,repl,cycle):
        """
        Returns true if an OpenMM replica has successfully completed a cycle.
        """
        try:
            pot = self._getPot(repl)
            if pot is None:
                return False
        except:
            return False
        return True

    def _getPar(self, repl):
        replica = self.openmm_replicas[repl]
        (stateid, par) = replica.get_state()
        return par

    def _getPot(self, repl):
        replica = self.openmm_replicas[repl]
        pot = replica.get_energy()
        return pot

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

        if self.keywords.get('LAMBDAS') is None:
            self._exit("LAMBDAS needs to be specified")
        self.lambdas = self.keywords.get('LAMBDAS').split(',')
        #list of temperatures
        if self.keywords.get('TEMPERATURES') is None:
            self._exit("TEMPERATURES needs to be specified")
        self.temperatures = self.keywords.get('TEMPERATURES').split(',')

        #flag to identify the intermediate states, typically the one at lambda=1/2
        self.intermediates = None
        self.intermediates = self.keywords.get('INTERMEDIATE').split(',')

        #direction of transformation at each lambda
        #ABFE 1 from RA to R+A, -1 from R+A to A
        #RBFE 1 from RA+B to RB+A, -1 from RB+A to RA+B
        self.directions = None
        self.directions = self.keywords.get('DIRECTION').split(',')

        #parameters of the softplus alchemical potential
        #lambda1 = lambda2 gives the linear potential
        self.lambda1s = None
        self.lambda2s = None
        self.alphas = None
        self.u0s = None
        self.w0coeffs = None
        self.lambda1s = self.keywords.get('LAMBDA1').split(',')
        self.lambda2s = self.keywords.get('LAMBDA2').split(',')
        self.alphas = self.keywords.get('ALPHA').split(',')
        self.u0s = self.keywords.get('U0').split(',')
        self.w0coeffs = self.keywords.get('W0COEFF').split(',')

        #build parameters for the lambda/temperatures combined states
        self.nreplicas = self._buildStates()

    def print_status(self):
        """
        Writes to BASENAME_stat.txt a text version of the status of the RE job

        It's fun to follow the progress in real time by doing:

        watch cat BASENAME_stat.txt
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
        direction = par['atmdirection']
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

    def _update_state_of_replica_addcustom(self, replica):
        #changes the format of the positions in case of an exchange between replicas with two different directions 
        #replica.convert_pos_into_direction_format()
        pass


class openmm_job_AmberRBFE(openmm_job_ATM):
    def __init__(self, command_file, options):
        super().__init__(command_file, options)

        prmtopfile = self.basename + ".prmtop"
        crdfile = self.basename + ".inpcrd"

        if self.stateparams is None:
            self._buildStates()

        #builds service worker for replicas use
        service_ommsys = OMMSystemAmberRBFE(self.basename, self.keywords, prmtopfile, crdfile, self.logger)
        self.service_worker = OMMWorkerATM(self.basename, service_ommsys, self.keywords, compute = False, logger = self.logger)
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
            ommsys = OMMSystemAmberRBFE(self.basename, self.keywords, prmtopfile, crdfile, self.logger) 
            self.openmm_workers.append(OMMWorkerATM(self.basename, ommsys, self.keywords, node_info = node, compute = True, logger = self.logger))
