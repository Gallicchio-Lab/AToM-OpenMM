import logging, logging.config
import math
import os

from configobj import ConfigObj
from openmm.unit import kelvin, kilocalories_per_mole

from gibbs_sampling import pairwise_independence_sampling
from ommreplica import OMMReplicaATM
from ommsystem import OMMSystemRBFE
from sync.worker import OMMWorkerATM
from utils.singal_guard import TerminationGuard
from utils.timer import Timer


class openmm_job_RBFE:

    def __init__(self, config_file):
        logging.config.fileConfig(os.path.join(os.path.dirname(__file__), "..", "utils", "logging.conf"))
        self.logger = logging.getLogger("sync_re")

        self.logger.info("Configuration:")
        self.config = ConfigObj(config_file)
        for key, value in self.config.items():
            self.logger.info(f"{key}: {value}")

        if self.config.get('VERBOSE').lower() == 'yes':
            self.logger.setLevel(logging.DEBUG)

        self.basename = self.config['BASENAME']
        self.state_params = self._getStateParams()
        self.nreplicas = len(self.state_params)

    def _getStateParams(self):
        lambdas = self.config['LAMBDAS'].split(',')
        directions = self.config['DIRECTION'].split(',')
        intermediates = self.config['INTERMEDIATE'].split(',')
        lambda1s = self.config['LAMBDA1'].split(',')
        lambda2s = self.config['LAMBDA2'].split(',')
        alphas = self.config['ALPHA'].split(',')
        u0s = self.config['U0'].split(',')
        w0s = self.config['W0COEFF'].split(',')
        temperatures = self.config['TEMPERATURES'].split(',')

        assert len(directions) == len(lambdas)
        assert len(intermediates) == len(lambdas)
        assert len(lambda1s) == len(lambdas)
        assert len(lambda2s) == len(lambdas)
        assert len(alphas) == len(lambdas)
        assert len(u0s) == len(lambdas)
        assert len(w0s) == len(lambdas)
        assert len(temperatures) == 1

        self.logger.info("State parameters")
        state_params = []
        for lambda_, direction, intermediate, lambda1, lambda2, alpha, u0, w0 in zip(lambdas, directions, intermediates, lambda1s, lambda2s, alphas, u0s, w0s):
            par = {}
            par['lambda'] = float(lambda_)
            par['atmdirection'] = float(direction)
            par['atmintermediate'] = float(intermediate)
            par['lambda1'] = float(lambda1)
            par['lambda2'] = float(lambda2)
            par['alpha'] = float(alpha)/kilocalories_per_mole
            par['u0'] = float(u0)*kilocalories_per_mole
            par['w0'] = float(w0)*kilocalories_per_mole
            par['temperature'] = float(temperatures[0])*kelvin
            state_params.append(par)
            self.logger.info(f"    State: {par}")

        return state_params

    def setupJob(self):
        with Timer(self.logger.info, "ATM setup"):

            with Timer(self.logger.info, "create system"):
                pdbtopfile = self.basename + ".pdb"
                systemfile = self.basename + "_sys.xml"
                ommsystem = OMMSystemRBFE(self.basename, self.config, pdbtopfile, systemfile, self.logger)
                ommsystem.create_system()

            with Timer(self.logger.info, "create worker"):
                self.worker = OMMWorkerATM(ommsystem, self.config, self.logger)

            with Timer(self.logger.info, "create replicas"):
                self.replicas = []
                for i in range(self.nreplicas):
                    replica = OMMReplicaATM(i, self.basename, self.worker, self.logger)
                    if not replica.get_stateid():
                        replica.set_state(i, self.state_params[i])
                    self.replicas.append(replica)

                self.replica_states = [replica.get_stateid() for replica in self.replicas]
                for i, replica in enumerate(self.replicas):
                    self.logger.info(f"Replica {i}: cycle {replica.get_cycle()}, state {replica.get_stateid()}")

            with Timer(self.logger.info, "update replicas"):
                self._updateReplicas()

    def scheduleJobs(self):
        with Timer(self.logger.info, "ATM simulations"):

            last_sample = self.replicas[0].get_cycle()
            num_samples = self.config['MAX_SAMPLES']
            if num_samples.startswith("+"):
                num_extra_samples = int(num_samples[1:])
                num_samples = num_extra_samples + last_sample - 1
                self.logger.info(f"Additional number of samples: {num_extra_samples}")
            else:
                num_samples = int(num_samples)
                self.logger.info(f"Target number of samples: {num_samples}")

            for isample in range(last_sample, num_samples + 1):
                with Timer(self.logger.info, f"sample {isample}"):

                    for irepl, replica in enumerate(self.replicas):
                        with Timer(self.logger.info, f"sample {isample}, replica {irepl}"):
                            assert replica.get_cycle() == isample
                            self.worker.run(replica)

                    with Timer(self.logger.info, "exchange replicas"):
                        self._exhangeReplicas()

                    with Timer(self.logger.info, "update replicas"):
                        self._updateReplicas()

                    with Timer(self.logger.info, "write replicas samples and trajectories"):
                        with TerminationGuard():
                            for replica in self.replicas:
                                replica.save_out()
                                if replica.get_mdsteps() % int(self.config['TRJ_FREQUENCY']) == 0:
                                    replica.save_dcd()

                    with Timer(self.logger.info, "checkpointing"):
                        with TerminationGuard():
                            for replica in self.replicas:
                                replica.save_checkpoint()

                    # Report progress on GPUGRID
                    progress = float(isample - last_sample + 1)/float(num_samples - last_sample + 1)
                    open("progress", "w").write(str(progress))

    def _updateReplicas(self):
        for replica, stateid in zip(self.replicas, self.replica_states):
            replica.set_state(stateid, self.state_params[stateid])

    def _exhangeReplicas(self):

        # Matrix of replica energies in each state.
        swap_matrix = self._computeSwapMatrix(range(self.nreplicas), self.replica_states)
        self.logger.debug("Swap matrix")
        for row in swap_matrix:
            self.logger.debug(f"    {row}")

        self.logger.debug(f"Replica states before: {self.replica_states}")
        for repl_i in range(self.nreplicas):
            sid_i = self.replica_states[repl_i]
            repl_j = pairwise_independence_sampling(repl_i,sid_i, range(self.nreplicas), self.replica_states, swap_matrix)
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
        U = [[ 0. for _ in range(self.nreplicas)] for _ in range(self.nreplicas)]

        n = len(repls)

        #collect replica parameters and potentials
        par = [self._getPar(k) for k in repls]
        pot = [self._getPot(k) for k in repls]

        for i in range(n):
            repl_i = repls[i]
            for j in range(n):
                sid_j = states[j]
                #energy of replica i in state j
                U[sid_j][repl_i] = self._reduced_energy(par[j],pot[i])
        return U

    def _getPar(self, repl):
        _, par = self.replicas[repl].get_state()
        return par

    #customized getPot to return the unperturbed potential energy
    #of the replica U0 = U - W_lambda(u)
    def _getPot(self, repl):
        replica = self.replicas[repl]
        pot = replica.get_energy()
        epot = pot['potential_energy']
        pertpot = pot['perturbation_energy']
        _, par = replica.get_state()
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

        kb = 0.0019872041*kilocalories_per_mole/kelvin
        beta = 1./(kb*par['temperature'])

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
