from __future__ import print_function
from __future__ import division
"""
Multiprocessing job transport for AsyncRE/OpenMM
"""
import os, re, sys, time, shutil, copy, random, signal
import multiprocessing as mp
#from multiprocessing import Process, Queue, Event
import logging

import openmm as mm
from openmm.app import *
from openmm import *
from openmm.unit import *
from datetime import datetime

from atom_openmm.ommreplica import *
from atom_openmm.ommsystem import *

from atom_openmm.utils.AtomUtils import AtomUtils

from contextlib import contextmanager

class OMMWorker(object):
    # OpenMM worker to run a replica in a process controlling one device
    #
    # The following methods can be overriden to support other replica exchange protocols.
    # See bedamtempt for an example
    #  set_state_values()
    #  get_energy_values()
    #  _worker_setstate_fromqueue()
    #  _worker_getenergy()
    #  _openmm_worker_body()
    def __init__(self, basename, ommsystem, keywords, node_info = None, compute = True, logger = None):
        self.node_name = None
        self.platform_name = None
        self.platformId = None
        self.deviceId = None
        self.nthreads = None
        if node_info is not None:
            self.node_name = node_info['node_name']
            pattern = re.compile(r'(\d+):(\d+)')
            self.platform_name = node_info['arch']
            matches = pattern.search(node_info["slot_number"])
            self.platformId = int(matches.group(1))
            self.deviceId = int(matches.group(2))
            self.nthreads = int(node_info["threads_number"])        
        self.basename = basename
        self.keywords = keywords
        self.ommsystem = ommsystem
        self.compute = compute
        self.logger = logger
        self.start_worker()

    def start_worker(self):
        self.ctx =  mp.get_context('spawn')
        self._startedSignal = self.ctx.Event()
        self._startedSignal.clear()
        self._readySignal = self.ctx.Event()
        self._readySignal.clear()
        self._runningSignal = self.ctx.Event()
        self._runningSignal.clear()
        self._errorSignal = self.ctx.Event()
        self._errorSignal.clear()
        self._isDone = self.ctx.Event()
        self._isDone.clear()

        self._cmdq = self.ctx.Queue()
        self._inq = self.ctx.Queue()
        self._outq = self.ctx.Queue()
        self.simulation = None
        self.context = None
        self.atmforce = None
        self.positions = None
        self.boxvectors = None
        self.par = {}
        self.pot = {}
        self.platform = None
        self.logfile_p = None
        self.outfile_p = None
        self.nprnt = int(self.keywords.get('PRNT_FREQUENCY'))
        if self.compute:
            #compute workers are launched as subprocesses
            s = signal.signal(signal.SIGINT, signal.SIG_IGN) #so that children do not respond to ctrl-c
            self._p = self.ctx.Process(target=self.openmm_worker, args=(self._startedSignal,self._readySignal,self._runningSignal,self._errorSignal,self._isDone,self._cmdq,self._inq,self._outq))
            self._p.daemon = True
            signal.signal(signal.SIGINT, s) #restore signal before start() of children
            self._p.start()
            self._readySignal.wait()
            return self._p
        else:
            #the service worker needs only the context in this process
            self._openmm_worker_body()
            self._openmm_worker_makecontext()
            return 1

    def set_state(self, par):
        self._readySignal.wait()
        self._cmdq.put("SETSTATE")
        self._inq.put(par)

    def get_energy(self):
        self._startedSignal.wait()
        self._readySignal.wait()
        self._cmdq.put("GETENERGY")
        pot = self._outq.get()
        return pot

    # set positions and velocities of worker
    def set_posvel(self, positions, velocities):
        self._startedSignal.wait()
        self._readySignal.wait()
        self._cmdq.put("SETPOSVEL")
        self._inq.put(positions)
        self._inq.put(velocities)

    # set positions and velocities of worker
    def set_chkpt(self, chkpt):
        self._startedSignal.wait()
        self._readySignal.wait()
        self._cmdq.put("SETCHKPT")
        self._inq.put(chkpt)
        
    # get positions and velocities from worker
    def get_posvel(self):
        self._startedSignal.wait()
        self._readySignal.wait()
        self._cmdq.put("GETPOSVEL")
        self.positions = self._outq.get()
        self.velocities = self._outq.get()
        return (self.positions, self.velocities)

    # set positions and velocities of worker
    def get_chkpt(self):
        self._startedSignal.wait()
        self._readySignal.wait()
        self._cmdq.put("GETCHKPT")
        chkpt = self._outq.get()
        return chkpt
    
    # sets the reporters of the worker
    def set_reporters(self, current_steps, outfile, logfile, xtcfile):
        self._startedSignal.wait()
        self._readySignal.wait()
        self._cmdq.put("SETREPORTERS")
        self._inq.put(current_steps)
        self._inq.put(logfile)

    # kills worker
    def finish(self, wait = False):
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

    # is worked done computing?
    def is_done(self):
        return self._isDone.is_set()

    # is worker started?
    def is_started(self):
        return self._startedSignal.is_set()

    # has worker died?
    def has_crashed(self):
        return not self._p.is_alive() or self._errorSignal.is_set()

    # starts execution loop of the worker
    def run(self, nsteps, nheating = 0, ncooling = 0, hightemp = 0.0):
        self._startedSignal.wait()
        self._readySignal.wait()
        self._cmdq.put("RUN")
        self._inq.put(nsteps)
        self._inq.put(nheating)
        self._inq.put(ncooling)
        self._inq.put(hightemp)
        self._runningSignal.set()
        self._isDone.clear()

    #
    # routine being multi-processed (worker)
    # default is temperature replica exchange (MD at constant temperature)
    #

    def _openmm_worker_body(self):
        
        self.ommsystem.create_system()
        self.system = self.ommsystem.system
        self.topology = self.ommsystem.topology
        self.integrator = self.ommsystem.integrator
        self.positions = self.ommsystem.positions
        self.boxvectors = self.ommsystem.boxvectors
        
    def _openmm_worker_run(self):
        try:
            if self.nheating > 0:
                self.integrator.setTemperature(self.hightemp)
                self.simulation.step(self.nheating)
                self.simulation.step(self.ncooling)
                production_temperature = self.par['temperature']
                self.integrator.setTemperature(production_temperature)
            self.simulation.step(self.nsteps)
            return 1
        except Exception as e:
            self.logger.error(f"MD has crashed: {e}")
            return None

    def _openmm_worker_makecontext(self):
        self.platform_properties = {}
        if self.platform_name is not None:
            if self.platform_name == 'OpenCL':
                self.platform = Platform.getPlatformByName(self.platform_name)
                self.platform_properties["OpenCLPlatformIndex"] = str(self.platformId)
                self.platform_properties["DeviceIndex"] = str(self.deviceId)
                self.platform_properties["Precision"] = "mixed"
                self.logger.info("Worker using OpenCL OpenMM platform")
            elif self.platform_name == "CUDA":
                self.platform = Platform.getPlatformByName(self.platform_name)
                self.platform_properties["DeviceIndex"] = str(self.deviceId)
                self.platform_properties["Precision"] = "mixed"
                self.logger.info("Worker using CUDA OpenMM platform")
            elif self.platform_name == "HIP":
                self.platform = Platform.getPlatformByName(self.platform_name)
                self.platform_properties["DeviceIndex"] = str(self.deviceId)
                self.platform_properties["Precision"] = "mixed"
                self.logger.info("Worker using HIP OpenMM platform")
            elif self.platform_name == "CPU":
                self.platform = Platform.getPlatformByName(self.platform_name)
                self.platform_properties["Threads"] = str(self.nthreads)
                self.logger.info("Worker using CPU OpenMM platform")
            elif self.platform_name == "Reference":
                self.platform = Platform.getPlatformByName(self.platform_name)
                self.logger.info("Worker using Reference OpenMM platform")
            else:
                self.logger.warning("Unrecognized platform name")
        else:
            self.platform = Platform.getPlatformByName('Reference')
            self.logger.info("Worker using Reference OpenMM platform")

        self.simulation = Simulation(self.topology, self.system, self.integrator, self.platform, self.platform_properties)
        self.context = self.simulation.context
        self.context.setPositions(self.positions)
        if self.boxvectors is not None:
            self.context.setPeriodicBoxVectors(*self.boxvectors)
        self.simulation.reporters = []

        #one preliminary energy evaluation seems to be required to init the energy routines
        if self.compute:
            self.simulation.context.applyConstraints(0.00001)
            state = self.simulation.context.getState(getEnergy = True)
            pote = state.getPotentialEnergy()

        #load initial state/coordinates
        self.simulation.loadState(self.basename + "_0.xml")

        #replace parameters loaded from the initial xml file with the values in the system
        for param_name in self.ommsystem.cparams:
            self.context.setParameter(param_name, self.ommsystem.cparams[param_name])
        
        if self.compute and self.platformId is not None and self.deviceId is not None:
            #sets up logfile
            self.wdir = "cntxt_%s_%d_%d" % (self.node_name,int(self.platformId),int(self.deviceId))
            if not os.path.isdir(self.wdir):
                os.mkdir(self.wdir)
            self.logfile = "%s/%s.log" % (self.wdir, self.basename)
            self.logfile_p = open(self.logfile, 'a+')
            self.simulation.reporters.append(StateDataReporter(self.logfile_p, self.nprnt, step=True, temperature=True, speed=True))

    def openmm_worker(self, startedSignal, readySignal, runningSignal, errorSignal, isDone, cmdq, inq, outq):
        try:
            import setproctitle
            setproctitle.setproctitle("AToM worker")
        except:
            pass
        
        startedSignal.clear()
        readySignal.clear()
        runningSignal.clear()
        errorSignal.clear()
        isDone.clear()

        self._openmm_worker_body()
        self._openmm_worker_makecontext()
        
        self.positions = None
        self.velocities = None
        self.chkpt = None

        #start event loop
        startedSignal.set()
        readySignal.set()
        while(True):
            readySignal.set()
            while(cmdq.empty()):
                time.sleep(0.1)
            command = cmdq.get()
            readySignal.clear()
            if command == "SETSTATE":
                self._worker_setstate(inq.get())
            elif command == "SETPOSVEL":
                self.positions = inq.get()
                self.velocities = inq.get()
                self.context.setPositions(self.positions)
                self.context.setVelocities(self.velocities)
            elif command == "RUN":
                runningSignal.set()
                isDone.clear()

                self.nsteps = int(inq.get())
                self.nheating = int(inq.get())
                self.ncooling = int(inq.get())
                self.hightemp = float(inq.get())

                #self.context.reinitialize(preserveState=True)
                
                res = self._openmm_worker_run()

                if self.logfile_p is not None:
                    self.logfile_p.flush()

                if res is None:
                    errorSignal.set()

                isDone.set()
            elif command == "GETENERGY":
                pot = self._worker_getenergy()
                outq.put(pot)
            elif command == "GETPOSVEL":
                state = self.context.getState(getPositions=True, getVelocities=True)
                self.positions = state.getPositions(asNumpy=True).value_in_unit(nanometers)
                self.velocities = state.getVelocities(asNumpy=True).value_in_unit(nanometers/picoseconds)
                outq.put(self.positions)
                outq.put(self.velocities)
            elif command == "SETCHKPT":
                self.chkpt = inq.get()
                self.context.loadCheckpoint(self.chkpt)
            elif command == "GETCHKPT":
                self.chkpt = self.context.createCheckpoint()
                outq.put(self.chkpt)
            elif command == "FINISH":
                if self.outfile_p is not None:
                    self.outfile_p.close()
                while not inq.empty():
                    inq.get()
                while not cmdq.empty():
                    cmdq.get()
                outq.close()
                cmdq.close()
                break
            else:
                #unknown command, do nothing, clear queues
                while not self.inq.empty():
                    inq.get()
                while not self.cmdq.empty():
                    cmdq.get()

        startedSignal.clear()
        readySignal.clear()

class OMMWorkerTRE(OMMWorker):
    def _worker_setstate(self, par):
        self.par = par
        self.integrator.setTemperature(self.par['temperature'])
        self.context.setParameter(self.ommsystem.parameter['temperature'], self.par['temperature']/kelvin)
  
    def _worker_getenergy(self):
        self.pot['potential_energy'] = self.context.getState(getEnergy = True).getPotentialEnergy()
        return self.pot

class OMMWorkerATM(OMMWorker):
    def _worker_setstate(self, par):
        self.par = par
        self.integrator.setTemperature(self.par['temperature'])
        self.context.setParameter(self.ommsystem.parameter['temperature'], self.par['temperature']/kelvin)
        atmforce = self.ommsystem.atmforce
        self.simulation.context.setParameter(atmforce.Lambda1(), self.par['lambda1'])
        self.simulation.context.setParameter(atmforce.Lambda2(), self.par['lambda2'])
        self.simulation.context.setParameter(atmforce.Alpha(), self.par['alpha']*kilojoules_per_mole)
        self.simulation.context.setParameter(atmforce.Uh(), self.par['uh'] /kilojoules_per_mole)
        self.simulation.context.setParameter(atmforce.W0(), self.par['w0'] /kilojoules_per_mole)
        self.simulation.context.setParameter(atmforce.Direction(), self.par['atmdirection'] )
        self.simulation.context.setParameter(atmforce.Umax(), self.par[atmforce.Umax()] /kilojoules_per_mole)
        self.simulation.context.setParameter(atmforce.Ubcore(), self.par[atmforce.Ubcore()] /kilojoules_per_mole)
        self.simulation.context.setParameter(atmforce.Acore(), self.par[atmforce.Acore()] )
        self.simulation.context.setParameter('UOffset', self.par['uoffset'] /kilojoules_per_mole )

    def _worker_getenergy(self):
        if self.ommsystem.v82plus:
            state = self.simulation.context.getState(getEnergy = True)
        else:
            if self.ommsystem.doMetaD:
                fgroups = { 0, self.ommsystem.metaDforcegroup, self.ommsystem.atmforcegroup }
            else:
                fgroups = { 0, self.ommsystem.atmforcegroup }
            state = self.simulation.context.getState(getEnergy = True, groups = fgroups )
        
        self.pot['potential_energy'] = state.getPotentialEnergy()
        (u1, u0, alchemicalEBias) = self.ommsystem.atmforce.getPerturbationEnergy(self.simulation.context)
        umcore = self.simulation.context.getParameter(self.ommsystem.atmforce.Umax())*kilojoules_per_mole
        ubcore = self.simulation.context.getParameter(self.ommsystem.atmforce.Ubcore())*kilojoules_per_mole
        acore = self.simulation.context.getParameter(self.ommsystem.atmforce.Acore())
        uoffset = self.simulation.context.getParameter('UOffset')*kilojoules_per_mole
        if self.par['atmdirection'] > 0:
            self.pot['perturbation_energy'] = self.ommsystem.atm_utils.softCorePertE(u1-(u0+uoffset), umcore, ubcore, acore)
        else:
            self.pot['perturbation_energy'] = self.ommsystem.atm_utils.softCorePertE((u0+uoffset)-u1, umcore, ubcore, acore)
        if self.ommsystem.doMetaD:
            state = self.simulation.context.getState(getEnergy = True, groups = {self.ommsystem.metaDforcegroup})
            self.pot['bias_energy'] = state.getPotentialEnergy()
        else:
            self.pot['bias_energy'] = 0.0 * kilojoules_per_mole
        return self.pot


class OMMWorkerATMSync(OMMWorkerATM):
    # Override start_worker to not add the queues and ctx
    def start_worker(self):
        self.simulation = None
        self.context = None
        self.atmforce = None
        self.positions = None
        self.velocities = None
        self.chkpt = None
        self.boxvectors = None
        self.par = {}
        self.pot = {}
        self.platform = None
        self.logfile_p = None
        self.outfile_p = None
        self.nprnt = int(self.keywords.get('PRNT_FREQUENCY'))
        #the service worker needs only the context in this process
        self._openmm_worker_body()
        self._openmm_worker_makecontext()
        return 1

    def set_state(self, par):
        self._worker_setstate(par)

    def get_energy(self):
        return self._worker_getenergy()

    # set positions and velocities of worker
    def set_posvel(self, positions, velocities):
        self.positions = positions
        self.velocities = velocities
        self.context.setPositions(self.positions)
        self.context.setVelocities(self.velocities)

    # set positions and velocities of worker
    def set_chkpt(self, chkpt):
        self.chkpt = chkpt
        self.context.loadCheckpoint(chkpt)
        
    # get positions and velocities from worker
    def get_posvel(self):
        state = self.context.getState(getPositions=True, getVelocities=True)
        self.positions = state.getPositions(asNumpy=True).value_in_unit(nanometers)
        self.velocities = state.getVelocities(asNumpy=True).value_in_unit(nanometers/picoseconds)
        return (self.positions, self.velocities)

    # set positions and velocities of worker
    def get_chkpt(self):
        self.chkpt = self.context.createCheckpoint()
        return self.chkpt
    
    # sets the reporters of the worker
    def set_reporters(self, current_steps, outfile, logfile, xtcfile):
        return

    # kills worker
    def finish(self, wait = False):
        return

    # is worker running?
    def is_running(self):
        return True

    # is worked done computing?
    def is_done(self):
        return True

    # is worker started?
    def is_started(self):
        return True

    # has worker died?
    def has_crashed(self):
        return False

    # starts execution loop of the worker
    def run(self, nsteps, nheating = 0, ncooling = 0, hightemp = 0.0):
        from atom_openmm.utils.timer import Timer

        self.nsteps = int(nsteps)
        self.nheating = int(nheating)
        self.ncooling = int(ncooling)
        self.hightemp = float(hightemp)

        #self.context.reinitialize(preserveState=True)
        with Timer(self.logger.info, "Executing replica"):
            self._openmm_worker_run()

        if self.logfile_p is not None:
            self.logfile_p.flush()
