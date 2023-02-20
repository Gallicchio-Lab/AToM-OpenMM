import logging
import math

from transport import Transport


class LocalOpenMMTransport(Transport):

    def __init__(self, jobname, worker, replicas):
        Transport.__init__(self)
        #self.logger = logging.getLogger("async_re.local_openmm_transport")
        self.logger = logging.getLogger("async_re.openmm_sync_re")

        self.worker = worker
        self.replicas = replicas

    def launchJob(self, replica, job_info):
        self.logger.debug('transport.lunchJob')

        _, par = replica.get_state()
        self.worker.set_state(par)
        self.worker.set_posvel(replica.positions, replica.velocities)

        self.worker.run(job_info['nsteps'])
        self._update_replica(replica, job_info)

    def _update_replica(self, replica, job_info):
        self.logger.debug(f'transport._update_replica: {job_info}')

        pos, vel = self.worker.get_posvel()
        assert pos
        assert vel

        pot = self.worker.get_energy()
        assert pot

        for value in pot.values():
            if math.isnan(value._value):
                return None
        for p in pos:
            if math.isnan(p.x) or math.isnan(p.y) or math.isnan(p.z):
                return None
        for v in vel:
            if math.isnan(v.x) or math.isnan(v.y) or math.isnan(v.z):
                return None

        replica.set_posvel(pos,vel)
        replica.set_energy(pot)

        cycle = replica.get_cycle() + 1
        replica.set_cycle(cycle)

        mdsteps = replica.get_mdsteps() + job_info['nsteps']
        replica.set_mdsteps(mdsteps)

        #output data and trajectory file update 
        if mdsteps % job_info['nprnt'] == 0:
            replica.save_out()
        if mdsteps % job_info['ntrj'] == 0:
            replica.save_dcd()
