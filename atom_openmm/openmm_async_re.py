from __future__ import print_function
from __future__ import division
import os
import re
import random
import math
import logging
import signal

import openmm as mm
from openmm.app import *
from openmm import *
from openmm.unit import *

from atom_openmm.async_re import JobManager
from atom_openmm.local_openmm_transport import *
from atom_openmm.ommreplica import *
from atom_openmm.ommsystem import *
from atom_openmm.ommworker import *

class openmm_job(JobManager):
    def __init__(self, command_file, options):
        super().__init__(command_file, options)
        self.openmm_replicas = None
        self.stateparams = None
        self.openmm_workers = None
        self.kb = 0.0019872041*kilocalories_per_mole/kelvin
        self.safeckpt_file = "ckpt_is_valid"
        
    def _setLogger(self):
        self.logger = logging.getLogger("async_re.openmm_async_re")
        
    def checkpointJob(self):
        #disable ctrl-c
        s = signal.signal(signal.SIGINT, signal.SIG_IGN)

        #the safeckpt_file is deleted at the start of the checkpointing
        #and created at the end. A missing checkpt_file when restarting
        #implies that the checkpointing stopped in the middle and that the
        #checkpoint files cannot be trusted
        try:
            os.remove(self.safeckpt_file)
        except:
            pass
            
        # update replica objects of waiting replicas
        self.update_replica_states()
        for replica in self.openmm_replicas:
            replica.save_checkpoint()

        #re-enable ctrl-c
        signal.signal(signal.SIGINT, s)

        #restore safeckpt_file
        with open(self.safeckpt_file, "w") as f:
            f.write("Checkpoint files are valid")

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
                replica.velocities *= scale

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

class openmm_job_TRE(openmm_job):
    def __init__(self, command_file, options):
        super().__init__(command_file, options)

        pdbtopfile = self.basename + ".pdb"
        systemfile = self.basename + "_sys.xml"

        if self.stateparams is None:
            self._buildStates()

        #builds service worker for replicas use
        service_ommsys = OMMSystemTRE(self.basename, self.keywords, pdbtopfile,  systemfile, self.logger)
        self.service_worker = OMMWorkerTRE(self.basename, ommsys, self.keywords, compute = False, logger = self.logger)
        #creates openmm replica objects
        self.openmm_replicas = []
        for i in range(self.nreplicas):
            try:
                replica = OMMReplicaTRE(i, self.basename, self.service_worker, self.logger, self.keywords)
            except:
                self._exit("Error creating replica.")
            if replica.stateid == None:
                replica.set_state(i, self.stateparams[i])#initial setting
            self.openmm_replicas.append(replica)

        # creates openmm workers
        self.openmm_workers = []
        pattern = re.compile(r'(\d+):(\d+)')
        for node in self.compute_nodes:
            slot_id = node["slot_number"]
            matches = pattern.search(slot_id)
            platform_id = int(matches.group(1))
            device_id = int(matches.group(2))
            gpu_platform_name = node["arch"]
            ommsys = OMMSystemTRE(self.basename, self.keywords,  pdbtopfile,  systemfile, self.logger)
            self.openmm_workers.append(OMMWorkerTRE(self.basename, ommsys, self.keywords, gpu_platform_name, platform_id, device_id, compute = True, logger = self.logger))

    def _buildStates(self):
        self.stateparams = []
        for tempt in self.temperatures:
            par = {}
            par['temperature'] = float(tempt)*kelvin
            self.stateparams.append(par)
        return len(self.stateparams)

    def _checkInput(self):
        super()._checkInput()

        if self.keywords.get('TEMPERATURES') is None:
            self._exit("TEMPERATURES needs to be specified")
        self.temperatures = self.keywords.get('TEMPERATURES')
        self.nreplicas = self._buildStates()

    def print_status(self):
        """
        Writes to BASENAME_stat.txt a text version of the status of the RE job

        It's fun to follow the progress in real time by doing:

        watch cat BASENAME_stat.txt
        """
        logfile = "%s_stat.txt" % self.basename
        ofile = open(logfile,"w")
        log = "Replica  State   Temperature Status  Cycle \n"
        for k in range(self.nreplicas):
            stateid = self.status[k]['stateid_current']
            log += "%6d   %5d  %s %6.2f  %5d \n" % (k, stateid, self.stateparams[stateid]['temperature']/kelvin, self.status[k]['running_status'], self.status[k]['cycle_current'])
        log += "Running = %d\n" % self.running
        log += "Waiting = %d\n" % self.waiting

        ofile.write(log)
        ofile.close()

    def _reduced_energy(self, par, pot):
        temperature = par['temperature']
        potential_energy = pot['potential_energy']
        beta = 1./(self.kb*temperature)
        return beta*potential_energy

        
class openmm_job_ATM(openmm_job):
    def _buildStates(self):
        self.stateparams = []
        for (lambd,direction,intermediate,lambda1,lambda2,alpha,uh,w0) in zip(self.lambdas,self.directions,self.intermediates,self.lambda1s,self.lambda2s,self.alphas,self.uhs,self.w0coeffs):
            for tempt in self.temperatures:
                par = {}
                par['lambda'] = float(lambd)
                par['atmdirection'] = float(direction)
                par['atmintermediate'] = float(intermediate)
                par['lambda1'] = float(lambda1)
                par['lambda2'] = float(lambda2)
                par['alpha'] = float(alpha)/kilocalories_per_mole
                par['uh'] = float(uh)*kilocalories_per_mole
                par['w0'] = float(w0)*kilocalories_per_mole
                par['temperature'] = float(tempt)*kelvin
                par['Umax'] = float(self.keywords.get('UMAX')) * kilocalorie_per_mole
                par['Ubcore'] = float(self.keywords.get('UBCORE')) * kilocalorie_per_mole
                par['Acore'] = float(self.keywords.get('ACORE'))
                if self.keywords.get('PERTE_OFFSET') is not None:
                    par['uoffset'] = float(self.keywords.get('PERTE_OFFSET')) * kilocalorie_per_mole
                else:
                    par['uoffset'] = 0.0 * kilocalorie_per_mole
                self.stateparams.append(par)
        return len(self.stateparams)

    def _checkInput(self):
        super()._checkInput()

        if self.keywords.get('LAMBDAS') is None:
            self._exit("LAMBDAS needs to be specified")
        self.lambdas = self.keywords.get('LAMBDAS')
        #list of temperatures
        if self.keywords.get('TEMPERATURES') is None:
            self._exit("TEMPERATURES needs to be specified")
        self.temperatures = self.keywords.get('TEMPERATURES')

        #flag to identify the intermediate states, typically the one at lambda=1/2
        self.intermediates = None
        self.intermediates = self.keywords.get('INTERMEDIATE')

        #direction of transformation at each lambda
        #ABFE 1 from RA to R+A, -1 from R+A to A
        #RBFE 1 from RA+B to RB+A, -1 from RB+A to RA+B
        self.directions = self.keywords.get('DIRECTION')

        #parameters of the softplus alchemical potential
        #lambda1 = lambda2 gives the linear potential
        self.lambda1s = None
        self.lambda2s = None
        self.alphas = None
        self.uhs = None
        self.w0coeffs = None
        self.lambda1s = self.keywords.get('LAMBDA1')
        self.lambda2s = self.keywords.get('LAMBDA2')
        self.alphas = self.keywords.get('ALPHA')
        self.uhs = self.keywords.get('U0')
        self.w0coeffs = self.keywords.get('W0COEFF')

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
        log = "Replica  State  Lambda Lambda1 Lambda2 Alpha Uh W0coeff Temperature Status  Cycle \n"
        for k in range(self.nreplicas):
            stateid = self.status[k]['stateid_current']
            log += "%6d   %5d  %6.3f %6.3f %6.3f %6.3f %6.2f %6.2f %6.2f %5s  %5d\n" % (k, stateid, self.stateparams[stateid]['lambda'], self.stateparams[stateid]['lambda1'], self.stateparams[stateid]['lambda2'], self.stateparams[stateid]['alpha']*kilocalories_per_mole, self.stateparams[stateid]['uh']/kilocalories_per_mole, self.stateparams[stateid]['w0']/kilocalories_per_mole, self.stateparams[stateid]['temperature']/kelvin, self.status[k]['running_status'], self.status[k]['cycle_current'])
        log += "Running = %d\n" % self.running
        log += "Waiting = %d\n" % self.waiting

        ofile.write(log)
        ofile.close()

    #evaluates the softplus function
    def _softplus(self, lambda1, lambda2, alpha, uh, w0, uf):
        ee = 1.0 + math.exp(-alpha*(uf-uh))
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
        uh = par['uh']
        w0 = par['w0']
        ebias = self._softplus(lambda1, lambda2, alpha, uh, w0, pertpot)
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
        uh = par['uh']
        w0 = par['w0']
        state_direction = par['atmdirection']
        state_intermediate = par['atmintermediate']
        epot0 = pot['unbiased_potential_energy']
        pertpot = pot['perturbation_energy']
        replica_direction = pot['direction']
        replica_intermediate = pot['intermediate']
        if (replica_direction == state_direction) or (state_intermediate > 0 and replica_intermediate > 0):
            ebias = self._softplus(lambda1, lambda2, alpha, uh, w0, pertpot)
            return beta*(epot0 + ebias)
        else:
            #prevent exchange
            large_energy = 1.e12
            return large_energy

    def _update_state_of_replica_addcustom(self, replica):
        #changes the format of the positions in case of an exchange between replicas with two different directions 
        #replica.convert_pos_into_direction_format()
        pass
            
class openmm_job_ABFE(openmm_job_ATM):
    def __init__(self, command_file, options):
        super().__init__(command_file, options)

        pdbtopfile = self.basename + ".pdb"
        systemfile = self.basename + "_sys.xml"

        if self.stateparams is None:
            self._buildStates()

        #builds service worker for replicas use
        service_ommsys = OMMSystemABFE(self.basename, self.keywords, pdbtopfile, systemfile, self.logger)
        self.service_worker = OMMWorkerATM(self.basename, service_ommsys, self.keywords, compute = False, logger = self.logger)
        #creates openmm replica objects
        self.openmm_replicas = []
        for i in range(self.nreplicas):
            try:
                replica = OMMReplicaATM(i, self.basename, self.service_worker, self.logger, self.keywords)
            except:
                self._exit("Error creating replica.")
            if replica.stateid == None:
                replica.set_state(i, self.stateparams[i])#initial setting
            self.openmm_replicas.append(replica)

        # creates openmm workers objects
        self.openmm_workers = []
        for node in self.compute_nodes:
            ommsys = OMMSystemABFE(self.basename, self.keywords, pdbtopfile, systemfile, self.logger) 
            self.openmm_workers.append(OMMWorkerATM(self.basename, ommsys, self.keywords, node_info = node, compute = True, logger = self.logger))

class openmm_job_RBFE(openmm_job_ATM):
    def __init__(self, command_file, options):
        super().__init__(command_file, options)

        self.async_mode = self.keywords.get('ASYNC_MODE', True)
        ommworkercls = OMMWorkerATM if self.async_mode else OMMWorkerATMSync
        
        pdbtopfile = self.basename + ".pdb"
        systemfile = self.basename + "_sys.xml"

        if self.stateparams is None:
            self._buildStates()

        #builds service worker for replicas use
        service_ommsys = OMMSystemRBFE(self.basename, self.keywords, pdbtopfile, systemfile, self.logger)
        self.service_worker = ommworkercls(self.basename, service_ommsys, self.keywords, compute = False, logger = self.logger)
        #creates openmm replica objects
        self.openmm_replicas = []
        for i in range(self.nreplicas):
            try:
                replica = OMMReplicaATM(i, self.basename, self.service_worker, self.logger, self.keywords)
            except:
                self._exit("Error when creating replica.")
            if replica.stateid == None:
                replica.set_state(i, self.stateparams[i])#initial setting
            self.openmm_replicas.append(replica)

        # creates openmm context objects
        self.openmm_workers = []
        for node in self.compute_nodes:
            ommsys = OMMSystemRBFE(self.basename, self.keywords, pdbtopfile, systemfile, self.logger) 
            self.openmm_workers.append(ommworkercls(self.basename, ommsys, self.keywords, node_info = node, compute = True, logger = self.logger))

