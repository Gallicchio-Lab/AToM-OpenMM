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

# routine being multi-processed
def openmm_worker(cmdq, outq, inq, startedSignal, readySignal, runningSignal,
                  basename, platform_id, device_id, keywords):
    startedSignal.clear()
    readySignal.clear()
    runningSignal.clear()

    rcptfile_input  = '%s_rcpt_0.dms' % basename 
    ligfile_input   = '%s_lig_0.dms'  % basename
        
    dms = DesmondDMSFile([ligfile_input, rcptfile_input]) 
    topology = dms.topology
    
    implicitsolvent = str(keywords.get('IMPLICITSOLVENT'))
    system = None
    if implicitsolvent is None:
        system = dms.createSystem(nonbondedMethod=NoCutoff, OPLS = True, implicitSolvent = None)
    elif implicitsolvent == 'AGBNP':
        system = dms.createSystem(nonbondedMethod=NoCutoff, OPLS = True, implicitSolvent = 'AGBNP')
    else:
        print('Unknown implicit solvent %s' % implicitsolvent)
        return

    natoms_ligand = int(keywords.get('NATOMS_LIGAND'))
    lig_atoms = range(natoms_ligand)
    # atom indexes here refer to indexes in either lig or rcpt dms file, rather than in the complex 
    #lig_atom_restr = [0, 1, 2, 3, 4, 5]   #indexes of ligand atoms for CM-CM Vsite restraint

    cm_lig_atoms = keywords.get('REST_LIGAND_CMLIG_ATOMS')   #indexes of ligand atoms for CM-CM Vsite restraint
    #convert the string of lig atoms to integer
    lig_atom_restr = [int(i) for i in cm_lig_atoms]
    #rcpt_atom_restr = [121, 210, 281, 325, 406, 527, 640, 650, 795, 976, 1276]   #indexes of rcpt atoms for CM-CM Vsite restraint
        
    cm_rcpt_atoms = keywords.get('REST_LIGAND_CMREC_ATOMS')   #indexes of rcpt atoms for CM-CM Vsite restraint
    #convert the string of receptor rcpt atoms to integer
    rcpt_atom_restr = [int(i) for i in cm_rcpt_atoms]

    cmkf = float(keywords.get('CM_KF'))
    kf = cmkf * kilocalorie_per_mole/angstrom**2 #force constant for Vsite CM-CM restraint
    cmtol = float(keywords.get('CM_TOL'))
    r0 = cmtol * angstrom #radius of Vsite sphere
    
    #these can be 'None" if not using orientational restraints
    lig_ref_atoms = None # the 3 atoms of the ligand that define the coordinate system of the ligand
    rcpt_ref_atoms = None # the 3 atoms of the receptor that define the coordinate system of the receptor
    angle_center = None * degrees
    kfangle = None * kilocalorie_per_mole/degrees**2
    angletol = None * degrees
    dihedral1center = None * degrees
    kfdihedral1 = None * kilocalorie_per_mole/degrees**2
    dihedral1tol = None * degrees
    dihedral2center = None * degrees
    kfdihedral2 = None * kilocalorie_per_mole/degrees**2
    dihedral2tol = None * degrees
    
    #transform indexes of receptor atoms
    for i in range(len(rcpt_atom_restr)):
        rcpt_atom_restr[i] += natoms_ligand
        if rcpt_ref_atoms:
            for i in range(len(rcpt_ref_atoms)):
                rcpt_ref_atoms[i] += natoms_ligand

    sdm_utils = SDMUtils(system, lig_atoms)
    sdm_utils.addRestraintForce(lig_cm_particles = lig_atom_restr,
                                rcpt_cm_particles = rcpt_atom_restr,
                                kfcm = kf,
                                tolcm = r0,
                                lig_ref_particles = lig_ref_atoms,
                                rcpt_ref_particles = rcpt_ref_atoms,
                                angle_center = angle_center,
                                kfangle = kfangle,
                                angletol = angletol,
                                dihedral1center = dihedral1center,
                                kfdihedral1 = kfdihedral1,
                                dihedral1tol = dihedral1tol,
                                dihedral2center = dihedral2center,
                                kfdihedral2 = kfdihedral2,
                                dihedral2tol = dihedral2tol)
    
    # the integrator object is context-specific
    #temperature = int(self.keywords.get('TEMPERATURES')) * kelvin
    temperature = 300 * kelvin
    frictionCoeff = float(keywords.get('FRICTION_COEFF')) / picosecond
    MDstepsize = float(keywords.get('TIME_STEP')) * picosecond
    umsc = float(keywords.get('UMAX')) * kilocalorie_per_mole
    acore = float(keywords.get('ACORE'))
    integrator = LangevinIntegratorSDM(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond, lig_atoms)
    integrator.setBiasMethod(sdm_utils.ILogisticMethod)
    integrator.setSoftCoreMethod(sdm_utils.RationalSoftCoreMethod)
    integrator.setUmax(umsc / kilojoule_per_mole)
    integrator.setAcore(acore)
    
    platform_properties = {}                                 
    platform = Platform.getPlatformByName('OpenCL')
    platform_properties["OpenCLPlatformIndex"] = str(platform_id)
    platform_properties["DeviceIndex"] = str(device_id)
    
    simulation = Simulation(topology, system, integrator, platform, platform_properties)
    context = simulation.context
    nsteps = int(keywords.get('PRODUCTION_STEPS'))
    nprnt = int(keywords.get('PRNT_FREQUENCY'))
    ntrj = int(keywords.get('TRJ_FREQUENCY'))
    if not nprnt % nsteps == 0:
        print("nprnt must be an integer multiple of nsteps.")
        return
    if not ntrj % nsteps == 0:
        print("ntrj must be an integer multiple of nsteps.")
        return
    simulation.reporters = []

    nsteps_current = 0
    outfile = None
    logfile = None
    dcdfile = None
    outfile_p = None
    logfile_p = None
    
    lmbd = None
    lmbd1 = None
    lmbd2 = None
    alpha = None
    u0 = None
    w0 = None

    pot_energy = None
    bind_energy = None
    
    positions = None
    velocities = None

    command = None

    write_out_flag = False
    write_trj_flag = False
    
    #start event loop
    startedSignal.set()
    readySignal.set()
    while(True):
        while(cmdq.empty()):
            time.sleep(0.1)
        command = cmdq.get()
        if command == "SETSTATE":
            lmbd = inq.get()
            lmbd1 = inq.get()
            lmbd2 = inq.get()
            alpha = inq.get()
            u0 = inq.get()
            w0 = inq.get()
            integrator.setLambda(lmbd)
            integrator.setLambda1(lmbd1)
            integrator.setLambda2(lmbd2)
            integrator.setAlpha(alpha*kilojoule_per_mole)
            integrator.setU0(u0/ kilojoule_per_mole)
            integrator.setW0coeff(w0 / kilojoule_per_mole)
        elif command == "SETPOSVEL":
            positions = inq.get()
            velocities = inq.get()
            context.setPositions(positions)
            context.setVelocities(velocities)
        elif command == "SETREPORTERS":
            replica_steps = inq.get()
            last_step = replica_steps + nsteps
            outfile = inq.get()
            logfile = inq.get()
            dcdfile = inq.get()
            simulation.reporters = []
            write_out_flag = ( last_step % nprnt == 0 )
            if last_step % nprnt == 0:
                if outfile_p != None:
                    outfile_p.close()
                outfile_p = open(outfile, 'a+')
                if logfile_p != None:
                    logfile_p.close()
                logfile_p = open(logfile, 'a+')
                simulation.reporters.append(StateDataReporter(logfile_p, nsteps, step=True, temperature=True))

            write_trj_flag = ( last_step % ntrj == 0 )
            if write_trj_flag :
                do_append = os.path.exists(dcdfile)
                simulation.reporters.append(DCDReporter(dcdfile, nsteps, append = do_append))
        elif command == "RUN":
            readySignal.clear()
            runningSignal.set()

            simulation.step(nsteps)

            nsteps_current += nsteps
            
            if write_out_flag:
                pot_energy = (integrator.getPotEnergy()*kilojoule_per_mole).value_in_unit(kilocalorie_per_mole)
                bind_energy = (integrator.getBindE()*kilojoule_per_mole).value_in_unit(kilocalorie_per_mole)
                if outfile_p:
                    outfile_p.write("%f %f %f %f %f %f %f %f\n" % (lmbd, lmbd1, lmbd2, alpha*kilocalorie_per_mole, u0/kilocalorie_per_mole, w0/kilocalorie_per_mole, pot_energy, bind_energy))
                    outfile_p.close()
                    outfile_p = None
                if logfile_p:
                    logfile_p.flush()

            simulation.reporters = [] #reset reporters
            
            runningSignal.clear()
            readySignal.set()
        elif command == "GETENERGY":
            bind_energy = (integrator.getBindE()*kilojoule_per_mole).value_in_unit(kilocalorie_per_mole)
            pot_energy = (integrator.getPotEnergy()*kilojoule_per_mole).value_in_unit(kilocalorie_per_mole)
            outq.put(bind_energy)
            outq.put(pot_energy)
        elif command == "GETPOSVEL":
            state = context.getState(getPositions=True, getVelocities=True)
            positions = state.getPositions()
            velocities = state.getVelocities()
            outq.put(positions)
            outq.put(velocities)
        elif command == "FINISH":
            startedSignal.clear()
            readySignal.clear()
            if outfile_p:
                outfile_p.close()
            break
        else:
            pass
        
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
        self._p = Process(target=openmm_worker, args=(self._cmdq,self._outq,self._inq, self._startedSignal, self._readySignal, self._runningSignal, basename, platform_id, device_id, keywords))
        signal.signal(signal.SIGINT, s) #restore signal before start() of children
        self._p.start()
    def set_state(self, lmbd, lmbd1, lmbd2, alpha, u0, w0):
        self.lmbd = float(lmbd)
        self.lmbd1 = float(lmbd1)
        self.lmbd2 = float(lmbd2)
        self.alpha = float(alpha)/kilocalorie_per_mole
        self.u0 = float(u0) * kilocalorie_per_mole
        self.w0 = float(w0) * kilocalorie_per_mole
        self._startedSignal.wait()
        self._readySignal.wait()
        self._cmdq.put("SETSTATE")
        self._inq.put(self.lmbd)
        self._inq.put(self.lmbd1)
        self._inq.put(self.lmbd2)
        self._inq.put(self.alpha)
        self._inq.put(self.u0)
        self._inq.put(self.w0)

    def get_energy(self):
        self._startedSignal.wait()
        self._readySignal.wait()
        self._cmdq.put("GETENERGY")
        bind_energy = self._outq.get()
        pot_energy = self._outq.get()
        return (bind_energy, pot_energy)
        
    def set_posvel(self, positions, velocities):
        self._startedSignal.wait()
        self._readySignal.wait()
        self._cmdq.put("SETPOSVEL")
        self._inq.put(positions)
        self._inq.put(velocities)

    def get_posvel(self):
        self._startedSignal.wait()
        self._readySignal.wait()
        self._cmdq.put("GETPOSVEL")
        self.positions = self._outq.get()
        self.velocities = self._outq.get()
        return (self.positions, self.velocities)

    def set_reporters(self, current_steps, outfile, logfile, dcdfile):
        self._startedSignal.wait()
        self._readySignal.wait()
        self._cmdq.put("SETREPORTERS")
        self._inq.put(current_steps)
        self._inq.put(outfile)
        self._inq.put(logfile)
        self._inq.put(dcdfile)
    
    def finish(self):
        self._cmdq.put("FINISH")
        self._p.join()

    def is_running(self):
        return self._runningSignal.is_set()

    def is_started(self):
        return self._startedSignal.is_set()
    
    def run(self):
        self._startedSignal.wait()
        self._readySignal.wait()
        self._cmdq.put("RUN")
        

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
        self.is_energy_assigned = False
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
        self.is_state_assigned = False
        
        # check for sdm_data table in lig dms
        tables = self.dms._tables[0]
        conn = self.sql_conn_lig
        if 'sdm_data' in tables:
            # read sdm_data table
            q = """SELECT binde,epot,temperature,lambda,lambda1,lambda2,alpha,u0,w0,cycle,stateid,mdsteps FROM sdm_data WHERE id = 1"""
            ans = conn.execute(q)
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
        self.is_state_assigned = True
        
    def get_state(self):
        return (self.bind_energy, self.pot_energy, self.stateid, self.lmbd, self.lambda1, self.lambda2, self.alpha, self.u0, self.w0coeff )

    def set_energy(self, bind_energy, pot_energy):
        self.bind_energy = bind_energy
        self.pot_energy = pot_energy
        self.is_energy_assigned = True
        
    def set_posvel(self, positions, velocities):
        self.positions = positions
        self.velocities = velocities

    def save_dms(self):
        if self.is_state_assigned and self.is_energy_assigned:
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
        sdm_context.set_state(replica.lmbd, replica.lambda1, replica.lambda2, replica.alpha, replica.u0, replica.w0coeff)
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
                (bind_energy, pot_energy ) = job['openmm_context'].get_energy()
                ommreplica.set_energy(bind_energy, pot_energy)
                #flag replica as not linked to a job
                self.replica_to_job[replica] = None

            return done
