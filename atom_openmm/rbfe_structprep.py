#! python

from __future__ import print_function
from __future__ import division
import sys
import time
import math
import random
import logging
import signal
import shutil
import random
from sys import stdout

import openmm as mm
from openmm.app import *
from openmm import *
from openmm.unit import *
from datetime import datetime

import logging

from atom_openmm.ommsystem import *
from atom_openmm.utils.AtomUtils import AtomUtils, residue_is_solvent

class OMMSystemRBFEnoATM(OMMSystemRBFE):
    def populate_atmforcegroup(self):
        import re
        nbpattern = re.compile(".*Nonbonded.*")
        gbpattern = re.compile(".*GB.*")
        harmpattern = re.compile(".*Harmonic.*")
        torpattern  = re.compile(".*Torsion.*")
        if self.var_force_group is not None:
            #add all forces in the var_force_group to ATMForce
            for i in range(self.system.getNumForces()):
                if self.system.getForce(i).getForceGroup() == self.var_force_group:
                    self.logger.info("Assigning Force %s to atmforcegroup" % self.system.getForce(i).getName())
                    self.system.getForce(i).setForceGroup(self.atmforcegroup)
        elif self.keywords.get('LIGAND1_VAR_ATOMS') is not None:
            #add all standard bonded and non-bonded forces to ATMForce
            for i in range(self.system.getNumForces()):
                if ( nbpattern.match(self.system.getForce(i).getName())   or
                     gbpattern.match(self.system.getForce(i).getName())   or
                     harmpattern.match(self.system.getForce(i).getName()) or
                     torpattern.match(self.system.getForce(i).getName())    ):
                    self.logger.info("Assigning Force %s to atmforcegroup" % self.system.getForce(i).getName())
                    self.system.getForce(i).setForceGroup(self.atmforcegroup)
        else:
            #add only non-bonded forces, after separating out the 1-4 interactions
            for i in range(self.system.getNumForces()):
                if ( nbpattern.match(self.system.getForce(i).getName()) or
                     gbpattern.match(self.system.getForce(i).getName()) ):
                    self.logger.info("Assigning Force %s to atmforcegroup" % self.system.getForce(i).getName())
                    self.system.getForce(i).setForceGroup(self.atmforcegroup)
                    
        self.logger.info("System's Forces:")
        for i in range(self.system.getNumForces()):
            self.logger.info("   %s" % self.system.getForce(i).getName())
                    
    def create_system(self):
        self.load_system()
        self.atm_utils = AtomUtils(self.system)
        self.set_ligand_atoms()
        self.set_displacement()
        self.set_vsite_restraints()
        #set orientation restraints
        self.set_orientation_restraints()
        #set reference atoms for alignment force
        self.set_alignmentForce()
        #indexes of the atoms whose position is restrained near the initial positions
        #by a flat-bottom harmonic potential.
        self.set_positional_restraints()
        
        #target temperature
        temp = self.keywords.get("TEMPERATURES")
        if temp == None:
            self.temperature = 300.0 * kelvin
        elif isinstance(temp, list):
            self.temperature = float(temp[0]) * kelvin
        else:
            self.temperature = float(temp) * kelvin

        self.set_torsion_metaDbias(self.temperature)

        #do not include ATM Force, instead place the nonbonded
        #forces in what it would be the the ATMForce group
        #for the integrator
        self.atmforcegroup = self.free_force_group()
        self.populate_atmforcegroup()

        #add barostat
        pressure=1*bar
        self.set_barostat(self.temperature,pressure,25)

        self.set_integrator(self.temperature, self.frictionCoeff, self.MDstepsize)

def set_platform(keywords):
    platform_properties = {}
    platform_name = keywords.get("OPENMM_PLATFORM")
    if platform_name == None:
        platform_name = "CUDA"
    if platform_name == "CUDA" or platform_name == "OpenCL" or platform_name == "HIP":
        platform_properties["Precision"] = "mixed"
    if platform_name == "CPU":
        try:
            nthreads = os.environ['OMP_NUM_THREADS']
        except:
            nthreads = 1
        platform_properties["Threads"] = str(nthreads)
    platform = Platform.getPlatformByName(platform_name)
    return (platform, platform_properties)

def do_mintherm(keywords, logger):
    basename = keywords.get('BASENAME')
    jobname = basename

    pdbtopfile = basename + ".pdb"
    systemfile = basename + "_sys.xml"

    #OpenMM system for minimization, thermalization, NPT, NVT
    #does not include ATM Force
    syst = OMMSystemRBFEnoATM(basename, keywords, pdbtopfile, systemfile, logger)
    syst.create_system()
    
    (platform, platform_properties) = set_platform(keywords)
    simulation = Simulation(syst.topology, syst.system, syst.integrator, platform, platform_properties)
    simulation.context.setPositions(syst.positions)
    if syst.boxvectors is not None:
        simulation.context.setPeriodicBoxVectors(syst.boxvectors[0], syst.boxvectors[1], syst.boxvectors[2] )
    simulation.context.applyConstraints(0.00001)
    print ("Using platform %s" % simulation.context.getPlatform().getName())
        
    print("Potential energy before minimization =", simulation.context.getState(getEnergy = True).getPotentialEnergy())
    print("Energy minimizing the system ...")
    simulation.minimizeEnergy()
    print("Potential energy after minimization =", simulation.context.getState(getEnergy = True).getPotentialEnergy())

    #saves minimization checkpoint
    simulation.saveState(jobname + '_min.xml')
    #saves a pdb file
    positions = simulation.context.getState(getPositions=True).getPositions()
    boxsize = simulation.context.getState().getPeriodicBoxVectors()
    simulation.topology.setPeriodicBoxVectors(boxsize)
    with open(jobname + '_min.pdb', 'w') as output:
        PDBFile.writeFile(simulation.topology, positions, output, keepIds=True)

    print("Thermalization ...")
 
    totalSteps = int(keywords.get("THERMALIZATION_STEPS", 150000))
    steps_per_cycle = int(keywords.get("STEPS_PER_CYCLE", 5000))
    number_of_cycles = int(totalSteps/steps_per_cycle)
    simulation.reporters.append(StateDataReporter(stdout, steps_per_cycle, step=True, potentialEnergy = True, temperature=True, volume=True, speed=True))

    # initial temperature
    initial_temp = keywords.get("INITIAL_TEMPERATURE")
    if initial_temp == None:
        initial_temperature = 50.0 * kelvin
    else:
        initial_temperature = float(initial_temp) * kelvin

    final_temperature = syst.temperature
    delta_temperature = (final_temperature - initial_temperature)/number_of_cycles

    syst.barostat.setFrequency(0)#disabled

    #MD with temperature ramp
    temperature = initial_temperature
    syst.integrator.setTemperature(temperature)
    for i in range(number_of_cycles):
        simulation.step(steps_per_cycle)
        #prepare system for new temperature
        temperature = temperature + delta_temperature
        syst.integrator.setTemperature(temperature)

    #saves thermalized checkpoint
    simulation.saveState(jobname + '_therm.xml')
    #saves a pdb file
    positions = simulation.context.getState(getPositions=True).getPositions()
    boxsize = simulation.context.getState().getPeriodicBoxVectors()
    simulation.topology.setPeriodicBoxVectors(boxsize)
    with open(jobname + '_therm.pdb', 'w') as output:
        PDBFile.writeFile(simulation.topology, positions, output, keepIds=True)

    print("NPT equilibration ...")
    
    syst.barostat.setFrequency(25)

    for i in range(number_of_cycles):
        simulation.step(steps_per_cycle)

    #saves checkpoint
    simulation.saveState(jobname + '_npt.xml')
    #saves a pdb file
    positions = simulation.context.getState(getPositions=True).getPositions()
    boxsize = simulation.context.getState().getPeriodicBoxVectors()
    simulation.topology.setPeriodicBoxVectors(boxsize)
    with open(jobname + '_npt.pdb', 'w') as output:
        PDBFile.writeFile(simulation.topology, positions, output, keepIds=True)

    print("NVT equilibration ...")
    syst.barostat.setFrequency(0)#disabled

    #MD at constant volume
    for i in range(number_of_cycles):
        simulation.step(steps_per_cycle)

    #saves checkpoint
    simulation.saveState(jobname + '_equil.xml')
    #saves a pdb file
    positions = simulation.context.getState(getPositions=True).getPositions()
    boxsize = simulation.context.getState().getPeriodicBoxVectors()
    simulation.topology.setPeriodicBoxVectors(boxsize)
    with open(jobname + '_equil.pdb', 'w') as output:
        PDBFile.writeFile(simulation.topology, positions, output, keepIds=True)

def do_lambda_annealing(keywords, logger):
    basename = keywords.get('BASENAME')
    jobname = basename
    
    pdbtopfile = basename + ".pdb"
    systemfile = basename + "_sys.xml"

    syst = OMMSystemRBFE(basename, keywords, pdbtopfile, systemfile, logger)
    syst.create_system()
    
    (platform, platform_properties) = set_platform(keywords)
    simulation = Simulation(syst.topology, syst.system, syst.integrator, platform, platform_properties)
    simulation.context.setPositions(syst.positions)
    if syst.boxvectors is not None:
        simulation.context.setPeriodicBoxVectors(syst.boxvectors[0], syst.boxvectors[1], syst.boxvectors[2])

    print ("Using platform %s" % simulation.context.getPlatform().getName())
    
    syst.barostat.setFrequency(0)#disabled

    #target temperature
    temp = keywords.get("TEMPERATURES")
    if temp == None:
        temperature = 300.0 * kelvin
    elif isinstance(temp, list):
        temperature = float(temp[0]) * kelvin
    else:
        temperature = float(temp) * kelvin

    lmbd = 0.0
    lambda1 = lmbd
    lambda2 = lmbd
    alpha = 0.0 / kilocalorie_per_mole
    uh = 0.0 * kilocalorie_per_mole
    w0coeff = 0.0 * kilocalorie_per_mole
    umsc =  1000.0 * kilocalorie_per_mole
    ubcore = 500.0 * kilocalorie_per_mole
    acore = 0.062500
    direction = 1
    uoffset = 0.0 * kilocalorie_per_mole

    print( "LoadState ...")
    simulation.loadState(jobname + '_equil.xml')

    #override ATM parameters
    simulation.context.setParameter(syst.atmforce.Lambda1(), lambda1)
    simulation.context.setParameter(syst.atmforce.Lambda2(), lambda2)
    simulation.context.setParameter(syst.atmforce.Alpha(), alpha *kilojoules_per_mole)
    simulation.context.setParameter(syst.atmforce.Uh(), uh /kilojoules_per_mole)
    simulation.context.setParameter(syst.atmforce.W0(), w0coeff /kilojoules_per_mole)
    simulation.context.setParameter(syst.atmforce.Umax(), umsc /kilojoules_per_mole)
    simulation.context.setParameter(syst.atmforce.Ubcore(), ubcore /kilojoules_per_mole)
    simulation.context.setParameter(syst.atmforce.Acore(), acore)
    simulation.context.setParameter(syst.atmforce.Direction(), direction)
    simulation.context.setParameter('UOffset',  uoffset /kilojoules_per_mole)

    epot = simulation.context.getState(getEnergy = True).getPotentialEnergy()
    print("Potential Energy =", simulation.context.getState(getEnergy = True, ).getPotentialEnergy())

    print("Annealing to lambda = 1/2 ...")

    totalSteps = int(keywords.get("ANNEALING_STEPS", 250000))
    steps_per_cycle = int(keywords.get("STEPS_PER_CYCLE", 5000))
    number_of_cycles = int(totalSteps/steps_per_cycle)
    deltalambda = (0.5 - 0.0)/float(number_of_cycles)
    simulation.reporters.append(StateDataReporter(stdout, steps_per_cycle, step=True, potentialEnergy = True, temperature=True, speed=True))
    if os.path.exists(jobname + "_mdlambda.xtc"):
        os.remove(jobname + "_mdlambda.xtc")
    simulation.reporters.append(XTCReporter(jobname + "_mdlambda.xtc", steps_per_cycle, enforcePeriodicBox=False))
    
    binding_file = jobname + '_mdlambda.out'
    f = open(binding_file, 'w')
    
    for i in range(number_of_cycles):
        simulation.step(steps_per_cycle)
        state = simulation.context.getState(getEnergy = True )
        pot_energy = state.getPotentialEnergy()
        (u1, u0, ebias) = syst.atmforce.getPerturbationEnergy(simulation.context)
        umcore = simulation.context.getParameter(syst.atmforce.Umax())* kilojoules_per_mole
        ubcore = simulation.context.getParameter(syst.atmforce.Ubcore())* kilojoules_per_mole
        acore = simulation.context.getParameter(syst.atmforce.Acore())
        direction = simulation.context.getParameter(syst.atmforce.Direction())
        if direction > 0:
            pert_energy = syst.atm_utils.softCorePertE(u1 - (u0+uoffset), umcore, ubcore, acore)
        else:
            pert_energy = syst.atm_utils.softCorePertE((u0+uoffset) - u1, umcore, ubcore, acore)
        l1 = simulation.context.getParameter(syst.atmforce.Lambda1())
        l2 = simulation.context.getParameter(syst.atmforce.Lambda2())
        a = simulation.context.getParameter(syst.atmforce.Alpha()) / kilojoules_per_mole
        umid = simulation.context.getParameter(syst.atmforce.Uh()) * kilojoules_per_mole
        w0 = simulation.context.getParameter(syst.atmforce.W0()) * kilojoules_per_mole
        print("%f %f %f %f %f %f %f %f %f" % (temperature/kelvin,lmbd, l1, l2, a*kilocalorie_per_mole, umid/kilocalorie_per_mole, w0/kilocalorie_per_mole, pot_energy/kilocalorie_per_mole, pert_energy/kilocalorie_per_mole), file=f )
        f.flush()
        lmbd += deltalambda
        lambda1 += deltalambda
        lambda2 += deltalambda
        simulation.context.setParameter(syst.atmforce.Lambda1(), lambda1)
        simulation.context.setParameter(syst.atmforce.Lambda2(), lambda2)

    f.close()
        
    print( "SaveState ...")
    simulation.saveState(jobname + "_mdlambda.xml")

    #saves a pdb file
    positions = simulation.context.getState(getPositions=True).getPositions()
    boxsize = simulation.context.getState().getPeriodicBoxVectors()
    simulation.topology.setPeriodicBoxVectors(boxsize)
    with open(jobname + '_mdlambda.pdb', 'w') as output:
        PDBFile.writeFile(simulation.topology, positions, output, keepIds=True)

def do_equil(keywords, logger):
    basename = keywords.get('BASENAME')
    jobname = basename

    pdbtopfile = basename + ".pdb"
    systemfile = basename + "_sys.xml"

    syst = OMMSystemRBFE(basename, keywords, pdbtopfile,  systemfile, logger)
    syst.create_system()

    (platform, platform_properties) = set_platform(keywords)
    simulation = Simulation(syst.topology, syst.system, syst.integrator, platform, platform_properties)
    simulation.context.setPositions(syst.positions)
    if syst.boxvectors is not None:
        simulation.context.setPeriodicBoxVectors(syst.boxvectors[0], syst.boxvectors[1], syst.boxvectors[2])

    print ("Using platform %s" % simulation.context.getPlatform().getName())
    
    syst.barostat.setFrequency(0)#disabled

    #target temperature
    temp = keywords.get("TEMPERATURES")
    if temp == None:
        temperature = 300.0 * kelvin
    elif isinstance(temp, list):
        temperature = float(temp[0]) * kelvin
    else:
        temperature = float(temp) * kelvin

    lmbd = 0.5
    lambda1 = lmbd
    lambda2 = lmbd
    alpha = 0.0 / kilocalorie_per_mole
    uh = 0.0 * kilocalorie_per_mole
    w0coeff = 0.0 * kilocalorie_per_mole
    umsc =  1000.0 * kilocalorie_per_mole
    ubcore = 500.0 * kilocalorie_per_mole
    acore = 0.062500
    direction = 1
    uoffset = 0.0 * kilocalorie_per_mole

    print( "LoadState ...")
    simulation.loadState(jobname + '_mdlambda.xml')

    #override ATM parameters
    simulation.context.setParameter(syst.atmforce.Lambda1(), lambda1)
    simulation.context.setParameter(syst.atmforce.Lambda2(), lambda2)
    simulation.context.setParameter(syst.atmforce.Alpha(), alpha *kilojoules_per_mole)
    simulation.context.setParameter(syst.atmforce.Uh(), uh /kilojoules_per_mole)
    simulation.context.setParameter(syst.atmforce.W0(), w0coeff /kilojoules_per_mole)
    simulation.context.setParameter(syst.atmforce.Umax(), umsc /kilojoules_per_mole)
    simulation.context.setParameter(syst.atmforce.Ubcore(), ubcore /kilojoules_per_mole)
    simulation.context.setParameter(syst.atmforce.Acore(), acore)
    simulation.context.setParameter(syst.atmforce.Direction(), direction)
    simulation.context.setParameter('UOffset',  uoffset /kilojoules_per_mole)

    epot = simulation.context.getState(getEnergy = True).getPotentialEnergy()
    print("Potential Energy =", simulation.context.getState(getEnergy = True).getPotentialEnergy())

    print("Equilibration at lambda = 1/2 ...")

    totalSteps = int(keywords.get("EQUILIBRATION_STEPS", 150000))
    steps_per_cycle = int(keywords.get("STEPS_PER_CYCLE", 5000))
    simulation.reporters.append(StateDataReporter(stdout, steps_per_cycle, step=True, potentialEnergy = True, temperature=True, speed=True))
    if os.path.exists(jobname + "_0.xtc"):
        os.remove(jobname + "_0.xtc")
    simulation.reporters.append(XTCReporter(jobname + "_0.xtc", steps_per_cycle, enforcePeriodicBox=False))
    
    simulation.step(totalSteps)
        
    print( "SaveState ...")
    simulation.saveState(jobname + "_0.xml")

    #saves a pdb file
    positions = simulation.context.getState(getPositions=True).getPositions()
    boxsize = simulation.context.getState().getPeriodicBoxVectors()
    simulation.topology.setPeriodicBoxVectors(boxsize)
    with open(jobname + '_0.pdb', 'w') as output:
        PDBFile.writeFile(simulation.topology, positions, output, keepIds=True)

def massage_keywords(keywords, restrain_solutes = True):

    #use 1 fs time step
    keywords['TIME_STEP'] = 0.001
    if keywords.get('MINTHERM_TIME_STEP') is not None:
        keywords['TIME_STEP'] = float(keywords.get('MINTHERM_TIME_STEP'))

    #temporarily restrain all non-solvent atoms
    if restrain_solutes:
        basename = keywords.get('BASENAME')
        pdbtopfile = basename + ".pdb"
        pdb = PDBFile(pdbtopfile)
        non_ion_wat_atoms = []
        for res in pdb.topology.residues():
            if not residue_is_solvent(res):
                for atom in res.atoms():
                    non_ion_wat_atoms.append(atom.index)
        keywords['POS_RESTRAINED_ATOMS'] = non_ion_wat_atoms


def rbfe_structprep(config_file=None):
    from atom_openmm.utils.AtomUtils import set_directory
    from atom_openmm.utils.config import parse_config
    from pathlib import Path

    if config_file is None:
        config_file = sys.argv[1]

    print("")
    print("========================================")
    print("AToM RBFE Structure Preparation         ")
    print("========================================")
    print("")
    print("Started at: " + str(time.asctime()))
    print("Input file:", config_file)
    print("")
    sys.stdout.flush()

    keywords = parse_config(config_file)
    logger = logging.getLogger("rbfe_structprep")

    with set_directory(Path(config_file).parent):
        restrain_solutes = True
        if keywords.get('MINTHERM_RESTRAIN_SOLUTES'):
            restrain_solutes = keywords.get('MINTHERM_RESTRAIN_SOLUTES').upper() == "YES"
        if not restrain_solutes:
            logger.info("Solutes are not restrained during mintherm.")
        old_keywords = keywords.copy()
        massage_keywords(keywords, restrain_solutes)

        do_mintherm(keywords, logger)
        do_lambda_annealing(keywords, logger)

        # reestablish the restrained atoms
        if restrain_solutes:
            keywords["POS_RESTRAINED_ATOMS"] = old_keywords.get("POS_RESTRAINED_ATOMS")

        do_equil(keywords, logger)


if __name__ == "__main__":
    assert len(sys.argv) == 2, "Specify ONE input file"

    rbfe_structprep(sys.argv[1])
