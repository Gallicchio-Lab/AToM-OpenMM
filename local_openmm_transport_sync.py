import logging
import math

from transport import Transport


class LocalOpenMMTransport(Transport):

    def __init__(self, jobname, openmm_workers, openmm_replicas):
        Transport.__init__(self)
        #self.logger = logging.getLogger("async_re.local_openmm_transport")
        self.logger = logging.getLogger("async_re.openmm_sync_re")

        # openmm contexts
        self.openmm_workers = openmm_workers
        assert len(self.openmm_workers) == 1

        # record replica OpenMM objects
        self.openmm_replicas = openmm_replicas

        # contains information about the device etc. running a replica
        # None = no information about where the replica is running
        self.replica_to_job = [ None for k in range(len(openmm_replicas)) ]

    def numNodesAlive(self):
        return 1

    def launchJob(self, replica, job_info):

        self.logger.debug('transport.lunchJob')
        job = job_info
        job['replica'] = replica
        self.replica_to_job[replica] = job

        node = 0
        job = self.replica_to_job[replica]

        job['nodeid'] = node
        job['openmm_replica'] = self.openmm_replicas[replica]
        job['openmm_worker'] = self.openmm_workers[node]

        self.replica_to_job[replica] = job
        self.LaunchReplica(job['openmm_worker'], job['openmm_replica'], job['cycle'],
                            job['nsteps'])

        nreplicas = len(self.replica_to_job)
        for repl in range(nreplicas):
            self.isDone(repl,0)

        return 1

    def LaunchReplica(self, worker, replica, cycle, nsteps):
        self.logger.debug('transport.LaunchReplica')
        _, par = replica.get_state()
        worker.set_posvel(replica.positions, replica.velocities)
        worker.set_state(par)
        worker.run(nsteps)

    def ProcessJobQueue(self, mintime, maxtime):
        return 1

    def DrainJobQueue(self):
        pass

    def _update_replica(self, job):
        self.logger.debug(f'transport._update_replica: {job}')

        #update replica cycle, mdsteps, write out, etc. from worker
        ommreplica = job['openmm_replica']
        pos, vel = job['openmm_worker'].get_posvel()
        assert pos
        assert vel

        pot = job['openmm_worker'].get_energy()
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

        cycle = ommreplica.get_cycle() + 1
        ommreplica.set_cycle(cycle)

        mdsteps = ommreplica.get_mdsteps() + job['nsteps']
        ommreplica.set_mdsteps(mdsteps)

        #update positions and velocities of openmm replica
        ommreplica.set_posvel(pos,vel)

        #update energies of openmm replica
        ommreplica.set_energy(pot)

        #output data and trajectory file update 
        if mdsteps % job['nprnt'] == 0:
            ommreplica.save_out()
        if mdsteps % job['ntrj'] == 0:
            ommreplica.save_dcd()
        return 0

    def isDone(self,replica,cycle):
        self.logger.debug(f'transport.isDone: {replica}')

        job = self.replica_to_job[replica]
        if job is None:
            # if job has been removed we assume that the replica is done
            return True

        if 'openmm_replica' not in job:
            return True

        self._update_replica(job)
        self.replica_to_job[replica] = None

        return True
