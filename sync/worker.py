import os
import sys

from openmm import OpenMMException, Platform
from openmm.app import Simulation, StateDataReporter
from openmm.unit import kelvin, kilojoules_per_mole

from utils.timer import Timer


class OMMWorkerATM:

    def __init__(self, ommsystem, config, logger):
        self.ommsystem = ommsystem
        self.config = config
        self.logger = logger

        self.topology = self.ommsystem.topology
        self.integrator = self.ommsystem.integrator

        nodefile = self.config.get('NODEFILE')
        assert nodefile, "NODEFILE needs to be specified"
        device = open(nodefile, 'r').readline().split(',')[1].strip().split(':')[1].strip()
        platform_name = open(nodefile, 'r').readline().split(',')[3].strip()

        if Platform.getNumPlatforms() == 1:
            if conda_prefix := os.environ.get("CONDA_PREFIX"):
                plugin_dir = os.path.join(conda_prefix, "lib", "plugins")
                Platform.loadPluginsFromDirectory(plugin_dir)

        platform = Platform.getPlatformByName(platform_name)
        properties = {"DeviceIndex": device, "Precision": "mixed"}
        self.logger.info(f"Device: {platform_name} {device}")

        self.simulation = Simulation(self.topology, self.ommsystem.system, self.integrator, platform, properties)
        self.context = self.simulation.context
        self.context.setPositions(self.ommsystem.positions)
        self.context.setPeriodicBoxVectors(*self.ommsystem.boxvectors)

        #one preliminary energy evaluation seems to be required to init the energy routines
        self.context.getState(getEnergy=True).getPotentialEnergy()

        #load initial state/coordinates
        basename = self.config['BASENAME']
        self.simulation.loadState(basename + "_0.xml")

        #replace parameters loaded from the initial xml file with the values in the system
        for key, value in self.ommsystem.cparams.items():
            self.context.setParameter(key, value)

        wdir = f"cntxt_{device}"
        if not os.path.isdir(wdir):
            os.mkdir(wdir)
        self.logfile = open(os.path.join(wdir, basename), 'a+')
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
        self.context.setParameter(atmforce.Uh(), par['uh'] /kilojoules_per_mole)
        self.context.setParameter(atmforce.W0(), par['w0'] /kilojoules_per_mole)
        self.context.setParameter(atmforce.Direction(), par['atmdirection'])

    def set_posvel(self, positions, velocities):
        self.logger.debug("ommworker.set_posvel")
        self.context.setPositions(positions)
        self.context.setVelocities(velocities)

    def get_energy(self, par):
        if self.ommsystem.doMetaD:
            fgroups = { 0, self.ommsystem.metaDforcegroup, self.ommsystem.atmforcegroup }
        else:
            fgroups = { 0, self.ommsystem.atmforcegroup }
        pot = {}
        state = self.context.getState(getEnergy = True, groups = fgroups)
        pot['potential_energy'] = state.getPotentialEnergy()
        
        (u1, u0, alchemicalEBias) = self.ommsystem.atmforce.getPerturbationEnergy(self.context)
        umcore = self.context.getParameter(self.ommsystem.atmforce.Umax())*kilojoules_per_mole
        ubcore = self.context.getParameter(self.ommsystem.atmforce.Ubcore())*kilojoules_per_mole
        acore = self.context.getParameter(self.ommsystem.atmforce.Acore())
        if par['atmdirection'] > 0:
            pot['perturbation_energy'] = self.ommsystem.atm_utils.softCorePertE(u1-u0, umcore, ubcore, acore)
        else:
            pot['perturbation_energy'] = self.ommsystem.atm_utils.softCorePertE(u0-u1, umcore, ubcore, acore)
        if self.ommsystem.doMetaD:
            state = self.simulation.context.getState(getEnergy = True, groups = {self.ommsystem.metaDforcegroup})
            pot['bias_energy'] = state.getPotentialEnergy()
        else:
            pot['bias_energy'] = 0.0 * kilojoules_per_mole
        self.logger.info(f'jack============> {pot}')
        return pot

    def get_posvel(self):
        self.logger.debug("ommworker.get_posvel")
        state = self.context.getState(getPositions=True, getVelocities=True)
        return state.getPositions(), state.getVelocities()

    def run(self, replica):
        assert replica.worker is self

        with Timer(self.logger.debug, "set replica state"):
            _, par = replica.get_state()
            self.set_state(par)
            self.set_posvel(replica.positions, replica.velocities)

        with Timer(self.logger.debug, "run replica"):
            nsteps = int(self.config['PRODUCTION_STEPS'])
            ntry = 5
            for _ in range(ntry):
                try:
                    self.simulation.step(nsteps)
                    break
                except OpenMMException as e:
                    self.logger.warning(f"Simulation failed: {e}")
            else:
                self.logger.error(f"Simulation failed {ntry} times!")
                sys.exit(0)

        with Timer(self.logger.debug, "get replica state"):
            pos, vel = self.get_posvel()
            _, par = replica.get_state()
            pot = self.get_energy(par)

            replica.set_posvel(pos, vel)
            replica.set_energy(pot)
            replica.set_cycle(replica.get_cycle() + 1)
            replica.set_mdsteps(replica.get_mdsteps() + nsteps)
