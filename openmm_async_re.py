from __future__ import print_function
from __future__ import division
import os
import re
import random
import math
import logging
import signal

from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *

from async_re import async_re
from local_openmm_transport import *
from ommreplica import *
from ommsystem import *
from ommworker import *

class openmm_job(async_re):
    def __init__(self, command_file, options):
        async_re.__init__(self, command_file, options)
        self.openmm_replicas = None
        self.stateparams = None
        self.openmm_workers = None
        self.gpu_platform_name = 'OpenCL' #default for now
        self.kb = 0.0019872041*kilocalories_per_mole/kelvin
        
    def _setLogger(self):
        self.logger = logging.getLogger("async_re.openmm_async_re")
        
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
        if not nprnt % nsteps == 0:
            self._exit("nprnt must be an integer multiple of nsteps.")
        if not ntrj % nsteps == 0:
            self._exit("ntrj must be an integer multiple of nsteps.")

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

    def _hasCompleted(self,replica,cycle):
        """
        Returns true if an OpenMM replica has successfully completed a cycle.
        """
        #safeguards are off for local transport
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
    def _buildStates(self):
        self.stateparams = []
        for tempt in self.temperatures:
            par = {}
            par['temperature'] = float(tempt)*kelvin
            self.stateparams.append(par)
        return len(self.stateparams)

    def _checkInput(self):
        async_re._checkInput(self)

        if self.keywords.get('TEMPERATURES') is None:
            self._exit("TEMPERATURES needs to be specified")
        self.temperatures = self.keywords.get('TEMPERATURES').split(',')
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
        for (lambd,lambda1,lambda2,alpha,u0,w0) in zip(self.lambdas,self.lambda1s,self.lambda2s,self.alphas,self.u0s,self.w0coeffs):
            for tempt in self.temperatures:
                par = {}
                par['lambda'] = float(lambd)
                par['lambda1'] = float(lambda1)
                par['lambda2'] = float(lambda2)
                par['alpha'] = float(alpha)/kilocalories_per_mole
                par['u0'] = float(u0)*kilocalories_per_mole
                par['w0'] = float(w0)*kilocalories_per_mole
                par['temperature'] = float(tempt)*kelvin
                self.stateparams.append(par)
        return len(self.stateparams)

    def _checkInput(self):
        async_re._checkInput(self)

        if self.keywords.get('LAMBDAS') is None:
            self._exit("LAMBDAS needs to be specified")
        self.lambdas = self.keywords.get('LAMBDAS').split(',')
        #list of temperatures
        if self.keywords.get('TEMPERATURES') is None:
            self._exit("TEMPERATURES needs to be specified")
        self.temperatures = self.keywords.get('TEMPERATURES').split(',')

        self.lambda1s = None
        self.lambda2s = None
        self.alphas = None
        self.u0s = None
        self.w0coeffs = None

        #ilogistic potential
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
        lambda1 = par['lambda1']
        lambda2 = par['lambda2']
        alpha = par['alpha']
        u0 = par['u0']
        w0 = par['w0']
        ebias = self._softplus(lambda1, lambda2, alpha, u0, w0, pertpot)
        pot['unbiased_potential_energy'] = epot - ebias
        return pot
        
    def _reduced_energy(self, par, pot):
        temperature = par['temperature']
        beta = 1./(self.kb*temperature)
        lambda1 = par['lambda1']
        lambda2 = par['lambda2']
        alpha = par['alpha']
        u0 = par['u0']
        w0 = par['w0']
        epot0 = pot['unbiased_potential_energy']
        pertpot = pot['perturbation_energy']
        ebias = self._softplus(lambda1, lambda2, alpha, u0, w0, pertpot)
        return beta*(epot0 + ebias)
        
class openmm_job_AmberTRE(openmm_job_TRE):
    def __init__(self, command_file, options):
        super().__init__(command_file, options)

        prmtopfile = self.basename + ".prmtop"
        crdfile = self.basename + ".inpcrd"

        if not self.stateparams:
            self._buildStates()

        #builds service worker for replicas use
        service_ommsys = OMMSystemAmberTRE(self.basename, self.keywords, prmtopfile, crdfile)
        self.service_worker = OMMWorkerTRE(self.basename, ommsys, self.keywords, compute = False)
        #creates openmm replica objects
        self.openmm_replicas = []
        for i in range(self.nreplicas):
            replica = OMMReplicaTRE(i, self.basename, self.service_worker)
            if replica.stateid == None:
                replica.set_state(i, self.stateparams[i])#initial setting
            self.openmm_replicas.append(replica)

        # creates openmm workers
        self.openmm_workers = []
        pattern = re.compile('(\d+):(\d+)')
        for node in self.compute_nodes:
            slot_id = node["slot_number"]
            matches = pattern.search(slot_id)
            platform_id = int(matches.group(1))
            device_id = int(matches.group(2))
            ommsys = OMMSystemAmberTRE(self.basename, self.keywords, prmtopfile, crdfile)
            self.openmm_workers.append(OMMWorkerTRE(self.basename, ommsys, self.keywords, self.gpu_platform_name, platform_id, device_id))
    
class openmm_job_AmberABFE(openmm_job_ATM):
    def __init__(self, command_file, options):
        super().__init__(command_file, options)

        prmtopfile = self.basename + ".prmtop"
        crdfile = self.basename + ".inpcrd"

        if not self.stateparams:
            self._buildStates()
        
        #builds service worker for replicas use
        service_ommsys = OMMSystemAmberABFE(self.basename, self.keywords, prmtopfile, crdfile)
        self.service_worker = OMMWorkerATM(self.basename, service_ommsys, self.keywords, compute = False)
        #creates openmm replica objects
        self.openmm_replicas = []
        for i in range(self.nreplicas):
            replica = OMMReplicaATM(i, self.basename, self.service_worker)
            if replica.stateid == None:
                replica.set_state(i, self.stateparams[i])#initial setting
            self.openmm_replicas.append(replica)

        # creates openmm workers objects
        self.openmm_workers = []
        pattern = re.compile('(\d+):(\d+)')
        for node in self.compute_nodes:
            slot_id = node["slot_number"]
            matches = pattern.search(slot_id)
            platform_id = int(matches.group(1))
            device_id = int(matches.group(2))
            ommsys = OMMSystemAmberABFE(self.basename, self.keywords, prmtopfile, crdfile) 
            self.openmm_workers.append(OMMWorkerATM(self.basename, ommsys, self.keywords, self.gpu_platform_name, platform_id, device_id))

class openmm_job_AmberRBFE(openmm_job_ATM):
    def __init__(self, command_file, options):
        super().__init__(command_file, options)

        prmtopfile = self.basename + ".prmtop"
        crdfile = self.basename + ".inpcrd"

        if not self.stateparams:
            self._buildStates()
        
        #builds service worker for replicas use
        service_ommsys = OMMSystemAmberRBFE(self.basename, self.keywords, prmtopfile, crdfile)
        self.service_worker = OMMWorkerATM(self.basename, service_ommsys, self.keywords, compute = False)
        #creates openmm replica objects
        self.openmm_replicas = []
        for i in range(self.nreplicas):
            replica = OMMReplicaATM(i, self.basename, self.service_worker)
            if replica.stateid == None:
                replica.set_state(i, self.stateparams[i])#initial setting
            self.openmm_replicas.append(replica)

        # creates openmm context objects
        self.openmm_workers = []
        pattern = re.compile('(\d+):(\d+)')
        for node in self.compute_nodes:
            slot_id = node["slot_number"]
            matches = pattern.search(slot_id)
            platform_id = int(matches.group(1))
            device_id = int(matches.group(2))
            ommsys = OMMSystemAmberRBFE(self.basename, self.keywords, prmtopfile, crdfile) 
            self.openmm_workers.append(OMMWorkerATM(self.basename, ommsys, self.keywords, self.gpu_platform_name, platform_id, device_id))



                                        
                                        
