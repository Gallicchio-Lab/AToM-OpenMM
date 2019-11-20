from __future__ import print_function
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

from transport import Transport

class OpenCLContext(object):
    # Context to run a SDM replica in a process controlling one OpenCL device
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
        return (pot_energy)
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
        self._inq.put(outfile)
        self._inq.put(logfile)
        self._inq.put(dcdfile)

    # kills worker
    def finish(self):
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

    # starts execution loop if the worker
    def run(self):
        self._startedSignal.wait()
        self._readySignal.wait()
        self._cmdq.put("RUN")

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
        pot_energy = (self.integrator.getPotEnergy()*kilojoule_per_mole).value_in_unit(kilocalorie_per_mole)
        temperature = self.par[0]
        if self.outfile_p:
            self.outfile_p.write("%f %f\n" % (temperature, pot_energy))

    def _worker_getenergy(self):                       
        pot_energy = (self.integrator.getPotEnergy()*kilojoule_per_mole).value_in_unit(kilocalorie_per_mole)
        self._outq.put(pot_energy)
        self.pot = [pot_energy]

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
        self.nsteps = int(self.keywords.get('PRODUCTION_STEPS'))
        self.nprnt = int(self.keywords.get('PRNT_FREQUENCY'))
        self.ntrj = int(self.keywords.get('TRJ_FREQUENCY'))
        if not self.nprnt % self.nsteps == 0:
            print("nprnt must be an integer multiple of nsteps.")
            return
        if not self.ntrj % self.nsteps == 0:
            print("ntrj must be an integer multiple of nsteps.")
            return
        self.simulation.reporters = []

        self.nsteps_current = 0
        self.outfile = None
        self.logfile = None
        self.dcdfile = None
        self.outfile_p = None
        self.logfile_p = None

        self.par = None
        self.pot = None
    
        self.positions = None
        self.velocities = None

        self.command = None

        write_out_flag = False
        write_trj_flag = False
    
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
            elif command == "SETREPORTERS":
                replica_steps = self._inq.get()
                last_step = replica_steps + self.nsteps
                outfile = self._inq.get()
                logfile = self._inq.get()
                dcdfile = self._inq.get()
                self.simulation.reporters = []
                write_out_flag = ( last_step % self.nprnt == 0 )
                if last_step % self.nprnt == 0:
                    if self.outfile_p != None:
                        self.outfile_p.close()
                    self.outfile_p = open(outfile, 'a+')
                    if self.logfile_p != None:
                        self.logfile_p.close()
                    self.logfile_p = open(logfile, 'a+')
                    self.simulation.reporters.append(StateDataReporter(self.logfile_p, self.nsteps, step=True, temperature=True))

                    write_trj_flag = ( last_step % self.ntrj == 0 )
                    if write_trj_flag :
                        do_append = os.path.exists(dcdfile)
                        self.simulation.reporters.append(DCDReporter(dcdfile, self.nsteps, append = do_append))
            elif command == "RUN":
                self._readySignal.clear()
                self._runningSignal.set()

                self.simulation.step(self.nsteps)

                self.nsteps_current += self.nsteps
            
                if write_out_flag:
                    self._worker_writeoutfile()
                    if self.outfile_p:
                        self.outfile_p.close()
                        self.outfile_p = None
                    if self.logfile_p:
                        self.logfile_p.flush()

                self.simulation.reporters = [] #reset reporters
            
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
        
class OMMReplica(object):
    #
    # Holds and manages OpenMM system for a replica
    #
    def __init__(self, replica_id, basename):
        self._id = replica_id
        self.basename = basename

        self.pot = None
        self.par = None
        self.is_energy_assigned = False
        self.is_state_assigned = False
        self.cycle = 0
        self.stateid = None
        self.mdsteps = 0
        
        self.open_dms()
        
        self.positions = copy.deepcopy(self.dms.positions)
        self.velocities = copy.deepcopy(self.dms.velocities)        


    def set_state(self, stateid, par):
        self.stateid = int(stateid)
        self.par = par
        self.is_state_assigned = True
        
    def get_state(self):
        return (self.stateid, self.par)

    def get_energy(self):
        return self.pot

    def set_energy(self, pot):
        self.pot = pot
        self.is_energy_assigned = True
        
    def set_posvel(self, positions, velocities):
        self.positions = positions
        self.velocities = velocities

    def open_dms(self):
        input_file  = '%s_0.dms' % self.basename 

        if not os.path.isdir('r%d' % self._id):
            os.mkdir('r%d' % self._id)

        output_file  = 'r%d/%s_ckp.dms' % (self._id,self.basename)
        if not os.path.isfile(output_file):
            shutil.copyfile(input_file, output_file)

        self.dms = DesmondDMSFile([output_file]) 

        self.sql_conn = self.dms._conn[0]
                
        # check for tre_data table in dms file
        tables = self.dms._tables[0]
        conn = self.sql_conn
        if 'tre_data' in tables:
            # read tre_data table
            q = """SELECT epot,temperature,cycle,stateid,mdsteps FROM tre_data WHERE id = 1"""
            ans = conn.execute(q)
            for (epot,temperature,cycle,stateid,mdsteps) in conn.execute(q):
                self.pot = (epot)
                self.par = (temperature)
                self.cycle = cycle
                self.stateid = stateid
                self.mdsteps = mdsteps
        else:
            #create tre_data table with dummy values
            conn.execute("CREATE TABLE IF NOT EXISTS tre_data (id INTEGER PRIMARY KEY, epot REAL, temperature REAL, cycle INTEGER, stateid INTEGER, mdsteps INTEGER )")
            conn.execute("INSERT INTO tre_data (epot,temperature,cycle,stateid,mdsteps) VALUES (0,0,0,0,0)")
            conn.commit()
            self.dms._tables[0] = self.dms._readSchemas(conn)

    def save_dms(self):
        if self.is_state_assigned and self.is_energy_assigned:
            conn = self.sql_conn
            pot_energy =  float(self.pot[0])            
            temperature = float(self.par[0])
            conn.execute("UPDATE tre_data SET epot = %f, temperature = %f, cycle = %d, stateid = %d, mdsteps = %d WHERE id = 1" % (pot_energy, temperature, self.cycle, self.stateid, self.mdsteps))
            conn.commit()
            self.dms.setPositions(self.positions)
            self.dms.setVelocities(self.velocities)

    def set_mdsteps(self, mdsteps):
        self.mdsteps = mdsteps

    def get_mdsteps(self):
        return self.mdsteps

    def set_cycle(self, cycle):
        self.cycle = cycle
        
    def get_cycle(self):
        return self.cycle

    def get_stateid(self):
        print(self.stateid)
        return self.stateid

class LocalOpenMMTransport(Transport):
    """
    Class to launch and monitor jobs on a set of local GPUs
    """
    def __init__(self, jobname, openmm_contexts, openmm_replicas):
        # jobname: identifies current asyncRE job
        # compute_nodes: list of OpenCL devices. The most important field is compute_nodes['slot_number']
        # which holds the device id in the format <platform_id>:<device_id>
        # nreplicas: number of replicas, 0 ... nreplicas-1
        Transport.__init__(self)
        self.logger = logging.getLogger("async_re.local_openmm_transport")
        
        # openmm contexts
        self.openmm_contexts = openmm_contexts
        self.nprocs = len(self.openmm_contexts)

        #constructs replica OpenMM objects
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

    def LaunchReplica(self, sdm_context, replica, cycle):
        (stateid, par) = replica.get_state()
        sdm_context.set_state(par)
        sdm_context.set_posvel(replica.positions, replica.velocities)
        
        #out_file     = 'r%d/%s_%d.out' % (replica._id,replica.basename,cycle)
        #log_file     = 'r%d/%s_%d.log' % (replica._id,replica.basename,cycle)
        #dcd_file     = 'r%d/%s_%d.dcd' % (replica._id,replica.basename,cycle)
        out_file     = 'r%d/%s.out' % (replica._id,replica.basename)
        log_file     = 'r%d/%s.log' % (replica._id,replica.basename)
        dcd_file     = 'r%d/%s.dcd' % (replica._id,replica.basename)
        sdm_context.set_reporters(replica.mdsteps,out_file, log_file, dcd_file)
        
        sdm_context.run()

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

                self.LaunchReplica(job['openmm_context'], job['openmm_replica'], job['cycle'])
                
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

            if not openmm_context.is_started():
                done = False
            else:
                done = not openmm_context.is_running()

            if done:
                # disconnects replica from job and node
                self._clear_resource(replica)
                #update cycle and mdsteps of replica
                ommreplica = job['openmm_replica']
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
                #flag replica as not linked to a job
                self.replica_to_job[replica] = None

            return done
