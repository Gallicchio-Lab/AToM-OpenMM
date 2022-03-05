from __future__ import print_function
from __future__ import division
"""
Multiprocessing job transport for AsyncRE/OpenMM
"""
import os, re, sys, time, shutil, copy, random, signal
import multiprocessing as mp
#from multiprocessing import Process, Queue, Event
import logging

from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from datetime import datetime

from ommreplica import *
from ommsystem import *

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
    def __init__(self, basename, ommsystem, keywords, platform_name = None , platform_id = None, device_id = None, compute = True, logger = None):
        self.platform_name = platform_name
        self.platformId = platform_id
        self.deviceId = device_id
        self.basename = basename
        self.keywords = keywords
        self.ommsystem = ommsystem
        self.compute = compute
        self.logger = logger
        self.start_worker()

    def start_worker(self):
        self.ctx =  mp.get_context('spawn')
        self._startedSignal = self.ctx.Event()
        self._readySignal = self.ctx.Event()
        self._runningSignal = self.ctx.Event()
        self._errorSignal = self.ctx.Event()
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
            self._p = self.ctx.Process(target=self.openmm_worker)
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
        except:
            self.logger.error("MD has crashed")
            return None

    def _openmm_worker_makecontext(self):
        self.platform_properties = {}
        if self.platform_name is not None:
            if self.platform_name == 'OpenCL':
                self.platform = Platform.getPlatformByName(self.platform_name)
                self.platform_properties["OpenCLPlatformIndex"] = str(self.platformId)
                self.platform_properties["DeviceIndex"] = str(self.deviceId)
            elif self.platform_name == "CUDA":
                self.platform = Platform.getPlatformByName(self.platform_name)
                self.platform_properties["DeviceIndex"] = str(self.deviceId)
            elif self.platform_name == "Reference":
                self.platform = Platform.getPlatformByName(self.platform_name)
            else:
                self.logger.warning("Unrecognized platform name")
        else:
            self.platform = Platform.getPlatformByName('Reference')

        self.simulation = Simulation(self.topology, self.system, self.integrator, self.platform, self.platform_properties)
        self.context = self.simulation.context
        self.context.setPositions(self.positions)
        if self.boxvectors is not None:
            self.context.setPeriodicBoxVectors(*self.boxvectors)
        self.simulation.reporters = []

        #one preliminary energy evaluation seems to be required to init the energy routines
        if self.compute:
            state = self.simulation.context.getState(getEnergy = True)#, groups = {1,3})
            pote = state.getPotentialEnergy()

        #load initial state/coordinates
        self.simulation.loadState(self.basename + "_0.xml")

        #replace parameters loaded from the initial xml file with the values in the system
        for param_name in self.ommsystem.cparams:
            self.context.setParameter(param_name, self.ommsystem.cparams[param_name])
        
        if self.compute and self.platformId is not None and self.deviceId is not None:
            #sets up logfile
            self.wdir = "cntxt%d_%d" % (int(self.platformId),int(self.deviceId))
            if not os.path.isdir(self.wdir):
                os.mkdir(self.wdir)
            self.logfile = "%s/%s.log" % (self.wdir, self.basename)
            self.logfile_p = open(self.logfile, 'a+')
            self.simulation.reporters.append(StateDataReporter(self.logfile_p, self.nprnt, step=True, temperature=True))

    def openmm_worker(self):
        
        self._startedSignal.clear()
        self._readySignal.clear()
        self._runningSignal.clear()
        self._errorSignal.clear()

        self._openmm_worker_body()
        self._openmm_worker_makecontext()
        
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

                res = self._openmm_worker_run()

                if self.logfile_p is not None:
                    self.logfile_p.flush()

                self._runningSignal.clear()
                self._readySignal.set()

                if res is None:
                    self._errorSignal.set()

            elif command == "GETENERGY":
                pot = self._worker_getenergy()
            elif command == "GETPOSVEL":
                state = self.context.getState(getPositions=True, getVelocities=True)
                self.positions = state.getPositions()
                self.velocities = state.getVelocities()
                self._outq.put(self.positions)
                self._outq.put(self.velocities)
            elif command == "FINISH":
                if self.outfile_p is not None:
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

        self._startedSignal.clear()
        self._readySignal.clear()

class OMMWorkerTRE(OMMWorker):
    def _worker_setstate_fromqueue(self):
        self.par = self._inq.get()
        self.integrator.setTemperature(self.par['temperature'])
        self.context.setParameter(self.ommsystem.parameter['temperature'], self.par['temperature']/kelvin)
  
    def _worker_getenergy(self):
        self.pot['potential_energy'] = self.context.getState(getEnergy = True).getPotentialEnergy()
        self._outq.put(self.pot)

class OMMWorkerATM(OMMWorker):
    def _worker_setstate_fromqueue(self):
        self.par = self._inq.get()
        self.integrator.setTemperature(self.par['temperature'])
        self.context.setParameter(self.ommsystem.parameter['temperature'], self.par['temperature']/kelvin)
        atmforce = self.ommsystem.atmforce
        self.simulation.context.setParameter(atmforce.Lambda1(), self.par['lambda1'])
        self.simulation.context.setParameter(atmforce.Lambda2(), self.par['lambda2'])
        self.simulation.context.setParameter(atmforce.Alpha(), self.par['alpha']*kilojoules_per_mole)
        self.simulation.context.setParameter(atmforce.U0(), self.par['u0'] /kilojoules_per_mole)
        self.simulation.context.setParameter(atmforce.W0(), self.par['w0'] /kilojoules_per_mole)
        self.simulation.context.setParameter(atmforce.Direction(), self.par['atmdirection'] )

    def _worker_getenergy(self):
        state = self.simulation.context.getState(getEnergy = True, groups = {1,3})
        self.pot['potential_energy'] = state.getPotentialEnergy()
        self.pot['perturbation_energy'] = self.ommsystem.atmforce.getPerturbationEnergy(self.simulation.context)
        self._outq.put(self.pot)
