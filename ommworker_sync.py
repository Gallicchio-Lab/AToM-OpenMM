import os
import re

from openmm import Platform
from openmm.app import Simulation, StateDataReporter
from openmm.unit import kelvin, kilojoules_per_mole


class OMMWorkerATM:

    def __init__(self, basename, ommsystem, keywords, node_info = None, compute = True, logger = None):
        self.node_name = None
        self.platform_name = None
        self.platformId = None
        self.deviceId = None
        self.nthreads = None
        if node_info is not None:
            self.node_name = node_info['node_name']
            pattern = re.compile('(\d+):(\d+)')
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

        self._openmm_worker_body()
        self._openmm_worker_makecontext()

    def set_state(self, par):
        self.logger.debug("ommworker.set_state")
        self.par = par
        self._worker_setstate_fromqueue()

    def get_energy(self):
        self.logger.debug("ommworker.get_energy")
        return self._worker_getenergy()

    def set_posvel(self, positions, velocities):
        self.logger.debug("ommworker.set_posvel")
        self.context.setPositions(positions)
        self.context.setVelocities(velocities)

    def get_posvel(self):
        self.logger.debug("ommworker.get_posvel")
        state = self.context.getState(getPositions=True, getVelocities=True)
        self.positions = state.getPositions()
        self.velocities = state.getVelocities()
        return (self.positions, self.velocities)

    def finish(self, wait = False):
        pass

    def is_running(self):
        return False

    def is_done(self):
        return True

    def is_started(self):
        return True

    def has_crashed(self):
        return False

    def run(self, nsteps, nheating = 0, ncooling = 0, hightemp = 0.0):
        self.logger.debug("ommworker.run")
        self.nsteps = nsteps
        self._openmm_worker_run()

    def _openmm_worker_body(self):
        self.ommsystem.create_system()
        self.system = self.ommsystem.system
        self.topology = self.ommsystem.topology
        self.integrator = self.ommsystem.integrator
        self.positions = self.ommsystem.positions
        self.boxvectors = self.ommsystem.boxvectors

    def _openmm_worker_run(self):
        self.simulation.step(self.nsteps)

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
            state = self.simulation.context.getState(getEnergy = True)#, groups = {1,3})
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
            self.simulation.reporters.append(StateDataReporter(self.logfile_p, self.nprnt, step=True, temperature=True))

    def _worker_setstate_fromqueue(self):

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

        fgroups = {0,self.ommsystem.metaDforcegroup,self.ommsystem.atmforcegroup}
        state = self.simulation.context.getState(getEnergy = True, groups = fgroups)

        self.pot['potential_energy'] = state.getPotentialEnergy()
        self.pot['perturbation_energy'] = self.ommsystem.atmforce.getPerturbationEnergy(self.simulation.context)
        self.pot['bias_energy'] = 0.0 * kilojoules_per_mole

        return self.pot
