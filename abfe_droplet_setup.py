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

from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from datetime import datetime

import logging
from configobj import ConfigObj
from atmmetaforce import *
from ommsystem import *


class OMMSystemAmberABFEnoATM(OMMSystemAmberABFE):
    def create_system(self):
        super().create_system(True, True, False)
        droplet = self.keywords.get('DROPLET').lower() == "true"
        if (droplet):
            print("Droplet calculation")
        else:
            msg = "Not a DROPLET calculation. ---- Use abfe_structprep.py"
            sys._exit(msg)
        self.load_amber_system(droplet)
        self.atm_utils = ATMMetaForceUtils(self.system)
        
        #target temperature
        temp = self.keywords.get("TEMPERATURES")
        if temp == None:
            self.temperature = 300.0 * kelvin
        elif isinstance(temp, list):
            self.temperature = float(temp[0]) * kelvin
        else:
            self.temperature = float(temp) * kelvin

        #self.set_torsion_metaDbias(self.temperature)

        #do not include ATM Force. 
        #self.set_atmforce()
        self.atmforcegroup = self.nonbondedforcegroup #for integrator




def do_mintherm(keywords, logger):
    basename = keywords.get('BASENAME')
    jobname = basename
    
    prmtopfile = basename + ".prmtop"
    crdfile = basename + ".inpcrd"

    #OpenMM system for minimization, thermalization, NVT
    #does not include ATM Force
    syst = OMMSystemAmberABFEnoATM(basename, keywords, prmtopfile, crdfile, logger)
    syst.create_system()
    
    platform_properties = {}
    platform_name = keywords.get("OPENMM_PLATFORM")
    if platform_name == None:
        platform_name = "CUDA"
    if platform_name == "CUDA" or platform_name == "OpenCL" or platform_name == "HIP":
        platform_properties["Precision"] = "mixed"
    platform = Platform.getPlatformByName(platform_name)
    
    simulation = Simulation(syst.topology, syst.system, syst.integrator, platform, platform_properties)
    simulation.context.setPositions(syst.positions)
    if syst.inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*syst.inpcrd.boxVectors)

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
        PDBFile.writeFile(simulation.topology, positions, output)

    print("Thermalization ...")

    #FIX ME - get from control file
    totalSteps = 100000
    steps_per_cycle = 5000
    number_of_cycles = int(totalSteps/steps_per_cycle)
    simulation.reporters.append(StateDataReporter(stdout, steps_per_cycle, step=True, potentialEnergy = True, temperature=True, volume=True))    
    simulation.reporters.append(DCDReporter(jobname + "_therm.dcd", steps_per_cycle))

    # initial temperature
    initial_temp = keywords.get("INITIAL_TEMPERATURE")
    if initial_temp == None:
        initial_temperature = 50.0 * kelvin
    else:
        initial_temperature = float(initial_temp) * kelvin

    final_temperature = syst.temperature
    delta_temperature = (final_temperature - initial_temperature)/number_of_cycles


    #MD with temperature ramp
    temperature = initial_temperature
    syst.integrator.setTemperature(temperature)
    for i in range(number_of_cycles):
        simulation.step(steps_per_cycle)
        #prepare system for new temperature
        temperature = temperature + delta_temperature
        syst.integrator.setTemperature(temperature)

    #saves thermalized checkpoint
    print( "SaveState ...")
    simulation.saveState(jobname + '_therm.xml')

    #saves a pdb file
    positions = simulation.context.getState(getPositions=True).getPositions()
    boxsize = simulation.context.getState().getPeriodicBoxVectors()
    simulation.topology.setPeriodicBoxVectors(boxsize)
    with open(jobname + '_therm.pdb', 'w') as output:
        PDBFile.writeFile(simulation.topology, positions, output)


def do_equil(keywords, logger):
    basename = keywords.get('BASENAME')
    jobname = basename
    
    prmtopfile = basename + ".prmtop"
    crdfile = basename + ".inpcrd"

    #syst = OMMSystemAmberABFE(basename, keywords, prmtopfile, crdfile, logger)
    syst = OMMSystemAmberABFEnoATM(basename, keywords, prmtopfile, crdfile, logger)
    syst.create_system()
    
    platform_properties = {}
    platform_name = keywords.get("OPENMM_PLATFORM")
    if platform_name == None:
        platform_name = "CUDA"
    if platform_name == "CUDA" or platform_name == "OpenCL" or platform_name == "HIP":
        platform_properties["Precision"] = "mixed"
    platform = Platform.getPlatformByName(platform_name)

    print("Equilibration ...")
    
    simulation = Simulation(syst.topology, syst.system, syst.integrator, platform, platform_properties)
    simulation.context.setPositions(syst.positions)
    if syst.inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*syst.inpcrd.boxVectors)

    print ("Using platform %s" % simulation.context.getPlatform().getName())

    print( "LoadState ...")
    simulation.loadState(jobname + '_therm.xml')

    if syst.doMetaD:
        fgroups = {0,syst.metaDforcegroup,syst.atmforcegroup}
    else:
        fgroups = {0,syst.atmforcegroup}


    state = simulation.context.getState(getEnergy = True, groups = fgroups)
    print("Potential Energy =", state.getPotentialEnergy())

    #FIX ME: get from keywords
    totalSteps = 100000
    steps_per_cycle = 5000
    simulation.reporters.append(StateDataReporter(stdout, steps_per_cycle, step=True, potentialEnergy = True, temperature=True, volume=True))
    simulation.reporters.append(DCDReporter(jobname + "_0.dcd", steps_per_cycle))

    state = simulation.context.getState(getEnergy = True, groups = fgroups)
    print("Potential Energy =", state.getPotentialEnergy())
    
    simulation.step(totalSteps)
        
    print( "SaveState ...")
    simulation.saveState(jobname + "_0.xml")

    #saves a pdb file
    positions = simulation.context.getState(getPositions=True).getPositions()
    boxsize = simulation.context.getState().getPeriodicBoxVectors()
    simulation.topology.setPeriodicBoxVectors(boxsize)
    with open(jobname + '_0.pdb', 'w') as output:
        PDBFile.writeFile(simulation.topology, positions, output)

def massage_keywords(keywords, restrain_solutes = False):

    #use 1 fs time step
    keywords['TIME_STEP'] = 0.001

    #restrain all solutes: receptor and ligands
    nlig = len(keywords.get('LIGAND_ATOMS'))
    last_lig_atom = int(keywords.get('LIGAND_ATOMS')[nlig-1])
    keywords['POS_RESTRAINED_ATOMS'] = [i for i in range(last_lig_atom+1)]

if __name__ == '__main__':

    # Parse arguments:
    usage = "%prog <ConfigFile>"

    if len(sys.argv) != 2:
        print("Please specify ONE input file")
        sys.exit(1)

    commandFile = sys.argv[1]

    print("")
    print("========================================")
    print("AToM ABFE Structure Preparation         ")
    print("========================================")
    print("")
    print("Started at: " + str(time.asctime()))
    print("Input file:", commandFile)
    print("")
    sys.stdout.flush()
    
    keywords = ConfigObj(commandFile)
    logger = logging.getLogger("rbfe_structprep")

    restrain_solutes = False
    old_keywords = keywords.copy()
    massage_keywords(keywords, restrain_solutes)
    
    do_mintherm(keywords, logger)

    #reestablish the restrained atoms
    if restrain_solutes:
        keywords['POS_RESTRAINED_ATOMS'] = old_keywords.get('POS_RESTRAINED_ATOMS') 

    do_equil(keywords, logger)

