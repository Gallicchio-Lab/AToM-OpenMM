import logging
import math

from transport import Transport


class LocalOpenMMTransport(Transport):

    def __init__(self, jobname, worker, config):
        Transport.__init__(self)
        #self.logger = logging.getLogger("async_re.local_openmm_transport")
        self.logger = logging.getLogger("async_re.openmm_sync_re")

        self.worker = worker
        self.config = config

    def launchJob(self, replica):
        self.logger.debug('transport.lunchJob')

        _, par = replica.get_state()
        self.worker.set_state(par)
        self.worker.set_posvel(replica.positions, replica.velocities)

        nsteps = int(self.config['PRODUCTION_STEPS'])
        self.worker.run(nsteps)

        pos, vel = self.worker.get_posvel()
        pot = self.worker.get_energy()

        replica.set_posvel(pos, vel)
        replica.set_energy(pot)
        replica.set_cycle(replica.get_cycle() + 1)
        replica.set_mdsteps(replica.get_mdsteps() + nsteps)
