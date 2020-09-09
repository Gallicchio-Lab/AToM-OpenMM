from __future__ import print_function
from __future__ import division
"""
Multiprocessing job transport for AsyncRE/OpenMM
"""
import os, re, sys, time, shutil, copy, random, signal
from multiprocessing import Process, Queue, Event
import logging

from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from simtk.openmm.app.desmonddmsfile import *
from datetime import datetime
from SDMplugin import *
from ommreplica import OMMReplica
from contextlib import contextmanager

from transport import Transport

class OpenCLContext(object):
    # Context to run a T-RE replica in a process controlling one OpenCL device
    #
    # The following methods can be overriden to support other replica exchange protocols.
    # See bedamtempt for an example
    #  set_state_values()
    #  get_energy_values()
    #  _worker_setstate_fromqueue()
    #  _worker_writeoutfile)
    #  _worker_getenergy()
    #  _openmm_worker_body()
    def __init__(self, basename, platform_id, device_id, keywords):
        s = signal.signal(signal.SIGINT, signal.SIG_IGN) #so that children do not respond to ctrl-c
        self._startedSignal = Event()
        self._readySignal = Event()
        self._runningSignal = Event()
        self._cmdq = Queue()
        self._inq = Queue()
        self._outq = Queue()
        self.platformId = platform_id
        self.deviceId = device_id
        self.basename = basename
        self.keywords = keywords
        self.nprnt = int(self.keywords.get('PRNT_FREQUENCY'))
        self._p = Process(target=self.openmm_worker)
        #, args=(self._cmdq,self._outq,self._inq, self._startedSignal, self._readySignal, self._runningSignal, basename, platform_id, device_id, keywords))
        self._p.daemon = True
        signal.signal(signal.SIGINT, s) #restore signal before start() of children
        self._p.start()

    # set the state of the worker
    def set_state_values(self, par):
        # can override in child class
        self.temperature = float(par[0])
        self._inq.put(self.temperature)
    def set_state(self, par):
        self._readySignal.wait()
        self._cmdq.put("SETSTATE")
        self.set_state_values(par)

    # get energies from the worker
    def get_energy_values(self):
        # can override in child class
        pot_energy = self._outq.get()
        return [pot_energy]
    def get_energy(self):
        self._startedSignal.wait()
        self._readySignal.wait()
        self._cmdq.put("GETENERGY")
        pot = self.get_energy_values()
        return pot

    # set positions and velocities of worker
    def set_posvel(self, positions, velocities):
        self._startedSignal.wait()
        self._readySignal.wait()
        self._cmdq.put("SETPOSVEL")
        self._inq.put(positions)
        self._inq.put(velocities)

    # get positions and velocities from worker
    def get_posvel(self):
        self._startedSignal.wait()
        self._readySignal.wait()
        self._cmdq.put("GETPOSVEL")
        self.positions = self._outq.get()
        self.velocities = self._outq.get()
        return (self.positions, self.velocities)

    # sets the reporters of the worker
    def set_reporters(self, current_steps, outfile, logfile, dcdfile):
        self._startedSignal.wait()
        self._readySignal.wait()
        self._cmdq.put("SETREPORTERS")
        self._inq.put(current_steps)
        self._inq.put(logfile)

    # kills worker
    def finish(self, wait = True):
        if wait:
            self._startedSignal.wait()
            self._readySignal.wait()
        self._cmdq.put("FINISH")
        while not self._outq.empty():
            self._outq.get()
        self._inq.close()
        self._cmdq.close()
        self._p.terminate()
        self._p.join(10) #10s time-out
        self._p.exitcode

    # is worker running?
    def is_running(self):
        return self._runningSignal.is_set()

    # is worker started?
    def is_started(self):
        return self._startedSignal.is_set()

    # has worker died?
    def has_crashed(self):
        return not self._p.is_alive()

    # starts execution loop of the worker
    def run(self, nsteps, nheating = 0, ncooling = 0, hightemp = 0.0):
        self._startedSignal.wait()
        self._readySignal.wait()
        self._cmdq.put("RUN")
        self._inq.put(nsteps)
        self._inq.put(nheating)
        self._inq.put(ncooling)
        self._inq.put(hightemp)

    #
    # routine being multi-processed (worker)
    # default is temperature replica exchange (MD at constant temperature)
    #
    def _worker_setstate_fromqueue(self):
        # can override in child class
        temperature = self._inq.get()
        self.integrator.setTemperature(temperature)
        self.par = [temperature]

    def _worker_writeoutfile(self):
        pot_energy = self.context.getState(getEnergy = True).getPotentialEnergy()/kilocalorie_per_mole
        temperature = self.par[0]
        if self.outfile_p:
            self.outfile_p.write("%f %f\n" % (temperature, pot_energy))

    def _worker_getenergy(self):
        pot_energy = self.context.getState(getEnergy = True).getPotentialEnergy()/kilocalorie_per_mole
        self._outq.put(pot_energy)
        self.pot = (pot_energy)

    def _openmm_worker_body(self):
        input_dms_file  = '%s_0.dms' % self.basename
        self.dms = DesmondDMSFile([input_dms_file])
        self.topology = self.dms.topology
        implicitsolvent = str(self.keywords.get('IMPLICITSOLVENT'))
        if implicitsolvent is None:
            self.system = self.dms.createSystem(nonbondedMethod=NoCutoff, OPLS = True, implicitSolvent = None)
        elif implicitsolvent == 'AGBNP':
            self.system = self.dms.createSystem(nonbondedMethod=NoCutoff, OPLS = True, implicitSolvent = 'AGBNP')
        else:
            print('Unknown implicit solvent %s' % implicitsolvent)
            sys.exit(1)

        temperature = 300 * kelvin #will be overridden by set_state()
        frictionCoeff = float(self.keywords.get('FRICTION_COEFF')) / picosecond
        MDstepsize = float(self.keywords.get('TIME_STEP')) * picosecond
        self.integrator = LangevinIntegrator(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond)

    def _openmm_worker_run(self):
        if self.nheating > 0:
            self.integrator.setTemperature(self.hightemp)
            self.simulation.step(self.nheating)
            self.simulation.step(self.ncooling)
            production_temperature = self.par[0]
            self.integrator.setTemperature(production_temperature)
        self.simulation.step(self.nsteps)

    def openmm_worker(self):
        self._startedSignal.clear()
        self._readySignal.clear()
        self._runningSignal.clear()

        self._openmm_worker_body()

        self.platform_properties = {}
        self.platform = Platform.getPlatformByName('OpenCL')
        self.platform_properties["OpenCLPlatformIndex"] = str(self.platformId)
        self.platform_properties["DeviceIndex"] = str(self.deviceId)

        self.simulation = Simulation(self.topology, self.system, self.integrator, self.platform, self.platform_properties)
        self.context = self.simulation.context

        self.simulation.reporters = []

        #sets up logfile
        self.wdir = "cntxt%d_%d" % (int(self.platformId),int(self.deviceId))
        if not os.path.isdir(self.wdir):
            os.mkdir(self.wdir)
        self.logfile = "%s/%s.log" % (self.wdir, self.basename)
        self.logfile_p = open(self.logfile, 'a+')
        self.simulation.reporters.append(StateDataReporter(self.logfile_p, self.nprnt, step=True, temperature=True))        

        self.par = None
        self.pot = None

        self.positions = None
        self.velocities = None

        self.command = None

        #start event loop
        self._startedSignal.set()
        self._readySignal.set()
        while(True):
            while(self._cmdq.empty()):
                time.sleep(0.1)
            command = self._cmdq.get()
            if command == "SETSTATE":
                self._worker_setstate_fromqueue()
            elif command == "SETPOSVEL":
                self.positions = self._inq.get()
                self.velocities = self._inq.get()
                self.context.setPositions(self.positions)
                self.context.setVelocities(self.velocities)
            elif command == "RUN":
                self._readySignal.clear()
                self._runningSignal.set()

                self.nsteps = int(self._inq.get())
                self.nheating = int(self._inq.get())
                self.ncooling = int(self._inq.get())
                self.hightemp = float(self._inq.get())

                self._openmm_worker_run()

                if self.logfile_p:
                    self.logfile_p.flush()

                self._runningSignal.clear()
                self._readySignal.set()
            elif command == "GETENERGY":
                pot = self._worker_getenergy()
            elif command == "GETPOSVEL":
                state = self.context.getState(getPositions=True, getVelocities=True)
                self.positions = state.getPositions()
                self.velocities = state.getVelocities()
                self._outq.put(self.positions)
                self._outq.put(self.velocities)
            elif command == "FINISH":
                if self.outfile_p:
                    self.outfile_p.close()
                while not self._inq.empty():
                    self._inq.get()
                while not self._cmdq.empty():
                    self._cmdq.get()
                self._outq.close()
                self._cmdq.close()
                break
            else:
                #unknown command, do nothing, clear queues
                while not self.inq.empty():
                    self._inq.get()
                while not self.cmdq.empty():
                    self._cmdq.get()

        startedSignal.clear()
        readySignal.clear()

class LocalOpenMMTransport(Transport):
    """
    Class to launch and monitor jobs on a set of local GPUs
    """
    def __init__(self, jobname, openmm_contexts, openmm_replicas):
        # jobname: identifies current asyncRE job
        Transport.__init__(self)
        self.logger = logging.getLogger("async_re.local_openmm_transport")

        # openmm contexts
        self.openmm_contexts = openmm_contexts
        self.nprocs = len(self.openmm_contexts)

        # record replica OpenMM objects
        self.openmm_replicas = openmm_replicas

        # device status = None if idle
        # Otherwise a structure containing:
        #    replica id being executed
        #    ...
        self.node_status = [ None for k in range(self.nprocs)]

        # contains information about the device etc. running a replica
        # None = no information about where the replica is running
        self.replica_to_job = [ None for k in range(len(openmm_replicas)) ]

        # implements a queue of jobs from which to draw the next job
        # to launch
        self.jobqueue = Queue()

    def _clear_resource(self, replica):
        # frees up the node running a replica identified by replica id
        job = None
        try:
            job = self.replica_to_job[replica]
        except:
            self.logger.warning("clear_resource(): unknown replica id %d",
                                replica)

        if job == None:
            return None

        try:
            nodeid = job['nodeid']
        except:
            self.logger.warning("clear_resource(): unable to query nodeid")
            return None

        try:
            if not self.node_status[nodeid] < 0: #signals a crashed node that should remain disabled
                self.node_status[nodeid] = None
        except:
            self.logger.warning("clear_resource(): unknown nodeid %", nodeid)
            return None

        return nodeid

    def _availableNode(self):
        #returns a node at random among available nodes
        available = [node for node in range(self.nprocs)
                     if self.node_status[node] == None]

        if available == None or len(available) == 0:
            return None
        random.shuffle(available)

        return available[0]

    def launchJob(self, replica, job_info):
        #Enqueues a replica for running based on provided job info.
        job = job_info
        job['replica'] = replica
        job['start_time'] = 0
        self.replica_to_job[replica] = job
        self.jobqueue.put(replica)
        return self.jobqueue.qsize()

    def LaunchReplica(self, sdm_context, replica, cycle, nsteps,
                      nheating = 0, ncooling = 0, hightemp = 0.0):
        (stateid, par) = replica.get_state()
        sdm_context.set_state(par)
        sdm_context.set_posvel(replica.positions, replica.velocities)
        sdm_context.run(nsteps, nheating, ncooling, hightemp)

    def ProcessJobQueue(self, mintime, maxtime):
        #Launches jobs waiting in the queue.
        #It will scan free devices and job queue up to maxtime.
        #If the queue becomes empty, it will still block until maxtime is elapsed.
        njobs_launched = 0
        usetime = 0
        nreplicas = len(self.replica_to_job)

        while usetime < maxtime:
            # find an available node
            node = self._availableNode()

            while (not self.jobqueue.empty()) and (not node == None):

                # grabs job on top of the queue
                replica = self.jobqueue.get()
                job = self.replica_to_job[replica]

                # assign job to available node
                job['nodeid'] = node
                job['openmm_replica'] = self.openmm_replicas[replica]
                job['openmm_context'] = self.openmm_contexts[node]
                job['start_time'] = time.time()

                # connects node to replica
                self.replica_to_job[replica] = job
                self.node_status[node] = replica

                if 'nheating' in job:
                    nheating = job['nheating']
                    ncooling = job['ncooling']
                    hightemp = job['hightemp']
                else:
                    nheating = 0
                    ncooling = 0
                    hightemp = 0.0

                self.LaunchReplica(job['openmm_context'], job['openmm_replica'], job['cycle'],
                                   job['nsteps'], nheating, ncooling, hightemp)

                # updates number of jobs launched
                njobs_launched += 1

                node = self._availableNode()

            # waits mintime second and rescans job queue
            time.sleep(mintime)

            # updates set of free nodes by checking for replicas that have exited
            for repl in range(nreplicas):
                self.isDone(repl,0)

            usetime += mintime

        return njobs_launched

    def DrainJobQueue(self):
        #clear the job queue
        while not self.jobqueue.empty():
            # grabs job on top of the queue
            replica = self.jobqueue.get()
            self._clear_resource(replica)
            self.replica_to_job[replica] = None

    def _update_replica(self, job):
        #update replica cycle, mdsteps, write out, etc. from context
        ommreplica = job['openmm_replica']
        if job['openmm_context'].has_crashed(): #refuses to update replica from a crashed context
            return
        cycle = ommreplica.get_cycle() + 1
        ommreplica.set_cycle(cycle)
        mdsteps = ommreplica.get_mdsteps() + job['nsteps']
        ommreplica.set_mdsteps(mdsteps)
        #update positions and velocities of openmm replica
        (pos,vel) = job['openmm_context'].get_posvel()
        ommreplica.set_posvel(pos,vel)
        #update energies of openmm replica
        pot = job['openmm_context'].get_energy()
        ommreplica.set_energy(pot)
        #output data and trajectory file update 
        if mdsteps % job['nprnt'] == 0:
            ommreplica.save_out()
        if mdsteps % job['ntrj'] == 0:
            ommreplica.save_dcd()

    def isDone(self,replica,cycle):
        """
        Checks if a replica completed a run.

        If a replica is done it clears the corresponding node.
        Note that cycle is ignored by job transport. It is assumed that it is
        the latest cycle.  it's kept for argument compatibility with
        hasCompleted() elsewhere.
        """
        job = self.replica_to_job[replica]
        if job == None:
            # if job has been removed we assume that the replica is done
            return True
        else:
            try:
                openmm_context = job['openmm_context']
            except:
                #job is in the queue but not yet launched
                return False

            if openmm_context.has_crashed():
                self.logger.warning("isDone(): replica %d has crashed", replica)
                openmm_context.finish(wait = False)
                self.node_status[job['nodeid']] = -1 #signals dead context
                self._clear_resource(replica)
                self.replica_to_job[replica] = None
                return True

            if not openmm_context.is_started():
                done = False
            else:
                done = not openmm_context.is_running()

            if done:
                # disconnects replica from job and node
                self._clear_resource(replica)
                #update replica info
                self._update_replica(job)
                #flag replica as not linked to a job
                self.replica_to_job[replica] = None

            return done
