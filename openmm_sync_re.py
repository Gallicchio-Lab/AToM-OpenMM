import logging, logging.config
import math
import os
import sys

from configobj import ConfigObj
from openmm.unit import kelvin, kilocalories_per_mole

from gibbs_sampling import pairwise_independence_sampling
from ommreplica import OMMReplicaATM
from ommsystem import OMMSystemAmberRBFE
from ommworker_sync import OMMWorkerATM
from utils.singal_guard import TerminationGuard
from utils.timer import Timer


class openmm_job_AmberRBFE:

    def __init__(self, config_file):
        logging.config.fileConfig(os.path.join(os.path.dirname(__file__), "utils/logging.conf"))
        self.logger = logging.getLogger("sync_re")

        self.logger.info("Configuration:")
        self.config = ConfigObj(config_file)
        for key, value in self.config.items():
            self.logger.info(f"{key}: {value}")

        if self.config.get('VERBOSE').lower() == 'yes':
            self.logger.setLevel(logging.DEBUG)

        self.basename = self.config['BASENAME']
        self.nreplicas = None
        self.kb = 0.0019872041*kilocalories_per_mole/kelvin

        self._checkInput()
        self._buildStates()

        # create system
        prmtopfile = self.basename + ".prmtop"
        crdfile = self.basename + ".inpcrd"
        ommsystem = OMMSystemAmberRBFE(self.basename, self.config, prmtopfile, crdfile, self.logger)
        ommsystem.create_system()

        # create worker
        self.worker = OMMWorkerATM(ommsystem, self.config, self.logger)

        #creates replicas
        self.openmm_replicas = []
        for i in range(self.nreplicas):
            replica = OMMReplicaATM(i, self.basename, self.worker, self.logger)
            if not replica.get_stateid():
                replica.set_state(i, self.stateparams[i])
            self.openmm_replicas.append(replica)

    def _checkInput(self):

        assert self.config.get('LAMBDAS'), "LAMBDAS needs to be specified"
        self.lambdas = self.config.get('LAMBDAS').split(',')
        assert self.config.get('TEMPERATURES'), "TEMPERATURES needs to be specified"
        self.temperatures = self.config.get('TEMPERATURES').split(',')
        assert len(self.temperatures) == 1

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

    def _buildStates(self):
        temperature = self.temperatures[0]
        self.stateparams = []
        for (lambd,direction,intermediate,lambda1,lambda2,alpha,u0,w0) in zip(self.lambdas,self.directions,self.intermediates,self.lambda1s,self.lambda2s,self.alphas,self.u0s,self.w0coeffs):
            par = {}
            par['lambda'] = float(lambd)
            par['atmdirection'] = float(direction)
            par['atmintermediate'] = float(intermediate)
            par['lambda1'] = float(lambda1)
            par['lambda2'] = float(lambda2)
            par['alpha'] = float(alpha)/kilocalories_per_mole
            par['u0'] = float(u0)*kilocalories_per_mole
            par['w0'] = float(w0)*kilocalories_per_mole
            par['temperature'] = float(temperature)*kelvin
            self.stateparams.append(par)
        return len(self.stateparams)

    def setupJob(self):
        # create status table
        self.replica_states = [replica.get_stateid() for replica in self.openmm_replicas]
        for i, replica in enumerate(self.openmm_replicas):
            self.logger.info(f"Replica {i}: cycle {replica.get_cycle()}, state {replica.get_stateid()}")

        self._updateReplicas()

    def scheduleJobs(self):

        assert 'MAX_SAMPLES' in self.config, "MAX_SAMPLES has to be specified"
        num_samples = int(self.config['MAX_SAMPLES'])

        last_sample = self.openmm_replicas[0].get_cycle()
        for isample in range(last_sample, num_samples + 1):
            with Timer(self.logger.info, f"sample {isample}"):

                for irepl, replica in enumerate(self.openmm_replicas):
                    with Timer(self.logger.info, f"sample {isample}, replica {irepl}"):
                        assert replica.get_cycle() == isample
                        self.worker.run(replica)

                with Timer(self.logger.info, "exchange replicas"):
                    self._exhangeReplicas()

                with Timer(self.logger.info, "update replicas"):
                    self._updateReplicas()

                with Timer(self.logger.info, "write replicas samples and trajectories"):
                    with TerminationGuard():
                        for replica in self.openmm_replicas:
                            replica.save_out()
                            replica.save_dcd()

                with Timer(self.logger.info, "checkpointing"):
                    with TerminationGuard():
                        for replica in self.openmm_replicas:
                            replica.save_checkpoint()

        self.logger.info("Done!")

    def _updateReplicas(self):
        for k in range(self.nreplicas):
            self._updateReplica(k)

    def _updateReplica(self, repl):
        replica = self.openmm_replicas[repl]
        old_stateid, _ = replica.get_state()
        new_stateid = self.replica_states[repl]
        replica.set_state(new_stateid, self.stateparams[new_stateid])

    def _exhangeReplicas(self):

        # Matrix of replica energies in each state.
        swap_matrix = self._computeSwapMatrix(range(self.nreplicas), self.replica_states)
        self.logger.debug("Swap matrix")
        for row in swap_matrix:
            self.logger.debug(f"    {row}")

        self.logger.debug(f"Replica states before: {self.replica_states}")
        for repl_i in range(self.nreplicas):
            sid_i = self.replica_states[repl_i]
            repl_j = pairwise_independence_sampling(repl_i,sid_i,
                                                    range(self.nreplicas),
                                                    self.replica_states,
                                                    swap_matrix)
            if repl_j != repl_i:
                sid_i = self.replica_states[repl_i]
                sid_j = self.replica_states[repl_j]
                self.replica_states[repl_i] = sid_j
                self.replica_states[repl_j] = sid_i
                self.logger.info(f"Replica {repl_i}: {sid_i} --> {sid_j}")
                self.logger.info(f"Replica {repl_j}: {sid_j} --> {sid_i}")

        self.logger.debug(f"Replica states after: {self.replica_states}")

    def _computeSwapMatrix(self, repls, states):
        """
        Compute matrix of dimension-less energies: each column is a replica
        and each row is a state so U[i][j] is the energy of replica j in state
        i.
        """
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

    def _getPar(self, repl):
        _, par = self.openmm_replicas[repl].get_state()
        return par

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
        # direction = par['atmdirection']
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

    #evaluates the softplus function
    def _softplus(self, lambda1, lambda2, alpha, u0, w0, uf):
        ee = 1.0 + math.exp(-alpha*(uf-u0))
        softplusf = lambda2 * uf + w0
        if alpha._value > 0.:
            softplusf += ((lambda2 - lambda1)/alpha) * math.log(ee)
        return softplusf
