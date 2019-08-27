from __future__ import print_function
"""
Multiprocessing job transport for AsyncRE/OpenMM
"""
import os, re, sys, time, shutil, copy, random
from threading import Thread, Event
import logging
import Queue

from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from simtk.openmm.app.desmonddmsfile import *
from datetime import datetime
from SDMplugin import *

from transport import Transport


class OpenCLContext(Thread):
    #
    # Context to run a SDM replica on one OpenCL device (reference platform also okay)
    #
    # based on reusable thread code from:
    # https://www.codeproject.com/Tips/1271787/Python-Reusable-Thread-Class
    def __init__(self, topology, system, integrator, platform,
                 platformId = None, deviceId = None , properties = None):
        self._startSignal = Event()
        self._oneRunFinished = Event()
        self._finishIndicator = False
        self.started = False
        self.running = False
        
        Thread.__init__(self)
        
        self.platform = platform

        if properties:
            self.platform_properties = properties
        else:
            self.platform_properties = {}
            
        if platformId:
            self.platformId = int(platformId)
            self.platform_properties["OpenCLPlatformIndex"] = str(platformId)
        else:
            self.platformId = None

        if deviceId:
            self.deviceId = int(deviceId)
            self.platform_properties["DeviceIndex"] = str(deviceId)
        else:
            self.deviceId = None


        #
        # I tried to make deep-copies of integrator, system,
        # etc. however it deepcopy() fails on the SDM and AGBNP plugin
        # modules complaining about "lack of proxy serialization". It
        # is probably due to some intricancy with swig. So, the
        # objects passed to this initializer should be created anew
        # for each opencl context, for example from the master .dms
        # files, otherwise they probably clobber each other.
        #
        self.topology = topology
        self.integrator = integrator
        self.system = system
        self.simulation = Simulation(self.topology, self.system, self.integrator, self.platform, self.platform_properties)
        self.context = self.simulation.context
        self.nsteps = None
        self.nprnt = None
        self.ntrj = None
        self.nsteps_current = 0
        self.outfile = None
        self.outfile_p = None
        
    def set_state(self, lmbd, lmbd1, lmbd2, alpha, u0, w0):
        self.lmbd = float(lmbd)
        self.lmbd1 = float(lmbd1)
        self.lmbd2 = float(lmbd2)
        self.alpha = float(alpha)/kilocalorie_per_mole
        self.u0 = float(u0) * kilocalorie_per_mole
        self.w0 = float(w0) * kilocalorie_per_mole
        self.integrator.setLambda(self.lmbd)
        self.integrator.setLambda1(self.lmbd1)
        self.integrator.setLambda2(self.lmbd2)
        self.integrator.setAlpha(self.alpha*kilojoule_per_mole)
        self.integrator.setU0(self.u0/ kilojoule_per_mole)
        self.integrator.setW0coeff(self.w0 / kilojoule_per_mole)

    def get_state(self):
        pot_energy = (self.integrator.getPotEnergy()*kilojoule_per_mole).value_in_unit(kilocalorie_per_mole)
        bind_energy = (self.integrator.getBindE()*kilojoule_per_mole).value_in_unit(kilocalorie_per_mole)
        return (bind_energy, pot_energy, self.lmbd, self.lmbd1, self.lmbd2, self.alpha, self.u0, self.w0 )
        
    def set_posvel(self, positions, velocities):
        self.context.setPositions(positions)
        self.context.setVelocities(velocities)

    def get_posvel(self):
        state = self.context.getState(getPositions=True, getVelocities=True)
        self.positions = state.getPositions()
        self.velocities = state.getVelocities()
        return (self.positions, self.velocities)

    def set_reporters(self, nprnt, outfile, logfile, ntrj, dcdfile):
        if not self.nsteps:
            self.nsteps = nprnt
        if not nprnt % self.nsteps == 0:
            print("nprnt must be an integer multiple of nsteps.")
            self.finish()
        else:
            self.nprnt = nprnt
        if not ntrj % self.nsteps == 0:
            print("ntrj must be an integer multiple of nsteps.")
            self.finish()
        else:
            self.ntrj = ntrj
        self.simulation.reporters = []
        last_step = self.nsteps_current + self.nsteps
        if last_step % self.nprnt == 0:
            self.simulation.reporters.append(StateDataReporter(logfile, nprnt, step=True, temperature=True))
            if self.outfile_p:
                self.outfile_p.close()
            self.outfile_p = open(outfile, 'w')
        if last_step % self.ntrj == 0:
            self.simulation.reporters.append(DCDReporter(dcdfile, ntrj))
        
    def outreport(self):
        (bind_energy, pot_energy, lmbd, lambda1, lambda2, alpha, u0, w0 ) = self.get_state()
        if self.outfile_p:
            self.outfile_p.write("%f %f %f %f %f %f %f %f\n" % (lmbd, lambda1, lambda2, alpha*kilocalorie_per_mole, u0/kilocalorie_per_mole, w0/kilocalorie_per_mole, pot_energy, bind_energy))
            self.outfile_p.flush()

    def set_nsteps(self, nsteps):
        if not self.nprnt:
            self.nprnt = nsteps
        else:
            if not self.nprnt % nsteps == 0:
                print("nprnt must be an integer multiple of nsteps.")
                self.finish()
        if not self.ntrj:
            self.ntrj = nsteps
        else:
            if not self.ntrj % nsteps == 0:
                print("ntrj must be an integer multiple of nsteps.")
                self.finish()
        self.nsteps = nsteps   
        
    def restart(self):
        """make sure to always call join() before restarting"""
        self._startSignal.set()

    def join(self):
        """ This join will only wait for one single run (target functioncall) to be finished"""
        self._oneRunFinished.wait()
        self._oneRunFinished.clear()

    def finish(self):
        self._finishIndicator = True
        self.restart()
        self.join()     

    def is_running(self):
        return self.running

    def is_started(self):
        return self.started
    
    def run(self):
        """ This class will reprocess the object "processObject" forever.
        Through the change of data inside processObject and start signals
        we can reuse the thread's resources"""

        self.started = True
        self.restart()
        while(True):    
            # wait until we should process
            self._startSignal.wait()

            self._startSignal.clear()

            if(self._finishIndicator):# check, if we want to stop
                self._oneRunFinished.set()
                return
            
            # call the threaded function
            self.running = True
            self.simulation.step(self.nsteps)
            self.nsteps_current += self.nsteps
            print("DEBUG nsteps_current: %d %d %d" % (self.nsteps, self.nsteps_current, self.nprnt))
            if self.nsteps_current % self.nprnt == 0:
                self.outreport()
            self.running = False

            # notify about the run's end
            self._oneRunFinished.set()

class SDMReplica(object):
    #
    # Holds and manages OpenMM system for a replica
    #
    def __init__(self, replica_id, basename):
        self._id = replica_id
        self.basename = basename
        
        rcptfile_input  = '%s_rcpt_0.dms' % basename 
        ligfile_input   = '%s_lig_0.dms'  % basename

        if not os.path.isdir('r%d' % replica_id):
            os.mkdir('r%d' % replica_id)

        ligfile_output  = 'r%d/%s_lig_ckp.dms' % (replica_id,basename)
        if not os.path.isfile(ligfile_output):
            shutil.copyfile(ligfile_input, ligfile_output)

        rcptfile_output = 'r%d/%s_rcpt_ckp.dms' % (replica_id,basename)
        if not os.path.isfile(rcptfile_output):
            shutil.copyfile(rcptfile_input, rcptfile_output)

        self.dms = DesmondDMSFile([ligfile_output, rcptfile_output]) 

        self.positions = copy.deepcopy(self.dms.positions)
        self.velocities = copy.deepcopy(self.dms.velocities)

        self.sql_conn_lig = self.dms._conn[0]
        self.sql_conn_rcpt = self.dms._conn[1]

        self.bind_energy = None
        self.pot_energy = None
        self.temperature = None
        self.lmbd =  None
        self.lambda1 =  None
        self.lambda2 =  None
        self.alpha = None
        self.u0 = None
        self.w0coeff = None
        self.cycle = 0
        self.stateid = None
        self.mdsteps = 0
        
        # check for sdm_data table in lig dms
        tables = self.dms._tables[0]
        conn = self.sql_conn_lig
        if 'sdm_data' in tables:
            # read sdm_data table
            print("reading checkpoint for replica %d" % replica_id)
            q = """SELECT binde,epot,temperature,lambda,lambda1,lambda2,alpha,u0,w0,cycle,stateid,mdsteps FROM sdm_data WHERE id = 1"""
            ans = conn.execute(q)
            print(ans)
            for (binde,epot,temperature,lmbd,lambda1,lambda2,alpha,u0,w0,cycle,stateid,mdsteps) in conn.execute(q):
                self.bind_energy = binde
                self.pot_energy = epot
                self.temperature = temperature
                self.lmbd =  lmbd
                self.lambda1 =  lambda1
                self.lambda2 =  lambda2
                self.alpha = alpha
                self.u0 = u0
                self.w0coeff = w0
                self.cycle = cycle
                self.stateid = stateid
                self.mdsteps = mdsteps
        else:
            #create sdm_data table with dummy values
            conn.execute("CREATE TABLE IF NOT EXISTS sdm_data (id INTEGER PRIMARY KEY, binde REAL, epot REAL, temperature REAL, lambda REAL, lambda1 REAL, lambda2 REAL, alpha REAL, u0 REAL, w0 REAL, cycle INTEGER, stateid INTEGER, mdsteps INTEGER )")
            conn.execute("INSERT INTO sdm_data (binde,epot,temperature,lambda,lambda1,lambda2,alpha,u0,w0,cycle,stateid,mdsteps) VALUES (0,0,0,0,0,0,0,0,0,0,0,0)")
            conn.commit()
            self.dms._tables[0] = self.dms._readSchemas(conn)

    def set_state(self, stateid, lmbd, lmbd1, lmbd2, alpha, u0, w0):
        self.stateid = int(stateid)
        self.lmbd = float(lmbd)
        self.lambda1 = float(lmbd1)
        self.lambda2 = float(lmbd2)
        self.alpha = float(alpha)
        self.u0 = float(u0)
        self.w0coeff = float(w0)
        self.temperature = 300. #not implemented

    def get_state(self):
        return (self.bind_energy, self.pot_energy, self.stateid, self.lmbd, self.lambda1, self.lambda2, self.alpha, self.u0, self.w0coeff )

    def set_energy(self, bind_energy, pot_energy):
        self.bind_energy = bind_energy
        self.pot_energy = pot_energy

    def set_posvel(self, positions, velocities):
        self.positions = positions
        self.velocities = velocities

    def save_dms(self):
        conn = self.sql_conn_lig
        conn.execute("UPDATE sdm_data SET binde = %f, epot = %f, temperature = %f, lambda = %f, lambda1 = %f, lambda2 = %f, alpha = %f, u0 = %f, w0 = %f, cycle = %d, stateid = %d, mdsteps = %d WHERE id = 1" % (self.bind_energy, self.pot_energy, self.temperature, self.lmbd, self.lambda1, self.lambda2, self.alpha, self.u0, self.w0coeff, self.cycle, self.stateid, self.mdsteps))
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
        self.jobqueue = Queue.Queue()

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

    def LaunchReplica(self, sdm_context, replica, cycle, nsteps, nprnt, ntrj):
        sdm_context.set_state(replica.lmbd, replica.lambda1, replica.lambda2, replica.alpha, replica.u0, replica.w0coeff)
        sdm_context.set_posvel(replica.positions, replica.velocities)
        
        out_file     = 'r%d/%s_%d.out' % (replica._id,replica.basename,cycle)
        log_file     = 'r%d/%s_%d.log' % (replica._id,replica.basename,cycle)
        dcd_file     = 'r%d/%s_%d.dcd' % (replica._id,replica.basename,cycle)
        sdm_context.set_nsteps(nsteps)
        sdm_context.set_reporters(nprnt, out_file, log_file, ntrj, dcd_file)
        
        if sdm_context.started:
            sdm_context.restart()
        else:
            sdm_context.start()

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

                print("DEBUG: Launching replica %d:" % replica)
                self.LaunchReplica(job['openmm_context'], job['openmm_replica'], job['cycle'], job['nsteps'], job['nprnt'], job['ntrj'])
                
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

    def isDone(self,replica,cycle):
        """
        Checks if a replica completed a run.

        If a replica is done it clears the corresponding node.
        Note that cycle is ignored by job transport. It is assumed that it is
        the latest cycle.  it's kept for argument compatibility with
        hasCompleted() elsewhere.
        """
        job = self.replica_to_job[replica]
        #print("DEBUG: isDone(): replica %d cycle %d" % (replica, cycle))
        #print(job)
        if job == None:
            # if job has been removed we assume that the replica is done
            return True
        else:
            try:
                openmm_context = job['openmm_context']
                print(job['openmm_context'])
            except:
                return False

            if  openmm_context == None:
                done = False
            else:
                if not openmm_context.is_started():
                    return False
                done = not openmm_context.is_running()

            if done:
                print("Replica %d completed cycle %d" % (replica, cycle))
                # disconnects replica from job and node
                self._clear_resource(replica)
                openmm_context.join()
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
                (bind_energy, pot_energy, lmbd, lambda1, lambda2, alpha, u0, w0 ) = job['openmm_context'].get_state()
                ommreplica.set_energy(bind_energy, pot_energy)
                #update checkpoint
                ommreplica.save_dms()
                #flag replica as not linked to a job
                self.replica_to_job[replica] = None

            return done
