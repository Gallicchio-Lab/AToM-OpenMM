import os
import re

from openmm import Platform
from openmm.app import Simulation, StateDataReporter
from openmm.unit import kelvin, kilojoules_per_mole


class OMMWorkerATM:

    def __init__(self, basename, ommsystem, config, logger=None):
        self.basename = basename
        self.ommsystem = ommsystem
        self.config = config
        self.logger = logger

        self.ommsystem.create_system()
        self.topology = self.ommsystem.topology
        self.integrator = self.ommsystem.integrator

        nodefile = self.config.get('NODEFILE')
        assert nodefile, "NODEFILE needs to be specified"
        device = open(nodefile, 'r').readline().split(',')[1].strip().split(':')[1].strip()

        platform = Platform.getPlatformByName("CUDA")
        properties = {"DeviceIndex": device, "Precision": "mixed"}
        self.logger.info(f"Device: CUDA {device}")

        self.simulation = Simulation(self.topology, self.ommsystem.system, self.integrator, platform, properties)
        self.context = self.simulation.context
        self.context.setPositions(self.ommsystem.positions)
        self.context.setPeriodicBoxVectors(*self.ommsystem.boxvectors)

        #one preliminary energy evaluation seems to be required to init the energy routines
        self.context.getState(getEnergy=True).getPotentialEnergy()

        #load initial state/coordinates
        self.simulation.loadState(self.basename + "_0.xml")

        #replace parameters loaded from the initial xml file with the values in the system
        for key, value in self.ommsystem.cparams.items():
            self.context.setParameter(key, value)

        self.wdir = f"cntxt_{device}"
        if not os.path.isdir(self.wdir):
            os.mkdir(self.wdir)
        self.logfile = open(os.path.join(self.wdir, self.basename), 'a+')
        nprnt = int(self.config.get('PRNT_FREQUENCY'))
        self.simulation.reporters.append(StateDataReporter(self.logfile, nprnt, step=True, temperature=True))

    def set_state(self, par):
        self.logger.debug("ommworker.set_state")

        self.integrator.setTemperature(par['temperature'])
        self.context.setParameter(self.ommsystem.parameter['temperature'], par['temperature']/kelvin)

        atmforce = self.ommsystem.atmforce
        self.context.setParameter(atmforce.Lambda1(), par['lambda1'])
        self.context.setParameter(atmforce.Lambda2(), par['lambda2'])
        self.context.setParameter(atmforce.Alpha(), par['alpha']*kilojoules_per_mole)
        self.context.setParameter(atmforce.U0(), par['u0'] /kilojoules_per_mole)
        self.context.setParameter(atmforce.W0(), par['w0'] /kilojoules_per_mole)
        self.context.setParameter(atmforce.Direction(), par['atmdirection'])

    def set_posvel(self, positions, velocities):
        self.logger.debug("ommworker.set_posvel")
        self.context.setPositions(positions)
        self.context.setVelocities(velocities)

    def run(self, nsteps):
        self.logger.info(f"Start MD simulation: {nsteps} steps")
        self.simulation.step(nsteps)
        self.logger.info("Finish MD simulation")

    def get_energy(self):
        self.logger.debug("ommworker.get_energy")
        fgroups = {0, self.ommsystem.atmforcegroup}
        state = self.context.getState(getEnergy = True, groups = fgroups)
        pot = {}
        pot['potential_energy'] = state.getPotentialEnergy()
        pot['perturbation_energy'] = self.ommsystem.atmforce.getPerturbationEnergy(self.context)
        pot['bias_energy'] = 0.0 * kilojoules_per_mole
        return pot

    def get_posvel(self):
        self.logger.debug("ommworker.get_posvel")
        state = self.context.getState(getPositions=True, getVelocities=True)
        return state.getPositions(), state.getVelocities()