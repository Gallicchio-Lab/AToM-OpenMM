from __future__ import print_function

from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os, re,time, shutil, math
from datetime import datetime

from atmmetaforce import *

print("Started at: " + str(time.asctime()))
start=datetime.now()

jobname = "<JOBNAME>"

lig_restr_offset = [  0.       for i in range(3) ] * angstrom


lig_atoms = [ <LIGATOMS> ]
rcpt_cm_atoms = [ <VSITERECEPTORATOMS> ]
restrained_atoms = [<RESTRAINEDATOMS>]

#load system
prmtop = AmberPrmtopFile(jobname + '.prmtop')
inpcrd = AmberInpcrdFile(jobname + '.inpcrd')
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                             constraints=HBonds)

#load the ATM Meta Force facility. Among other things the initializer
#sorts Forces into groups 
atm_utils = ATMMetaForceUtils(system)


lig_cm_atoms = lig_atoms
#vsite_restraints
kf = 25.0 * kilocalorie_per_mole/angstrom**2
r0 = 5.0 * angstrom # change it to s placeholder

atm_utils.addRestraintForce(lig_cm_particles = lig_cm_atoms,
                            rcpt_cm_particles = rcpt_cm_atoms,
                            kfcm = kf,
                            tolcm = r0,
                            offset = lig_restr_offset)

#receptor positional restraints
fc = 25.0 * kilocalorie_per_mole/angstrom**2
tol = 1.5 * angstrom

atm_utils.addPosRestraints(restrained_atoms, inpcrd.positions, fc, tol)

#define a barostat but disable it
#Set up Langevin integrator
temperature = 300 * kelvin
frictionCoeff = 0.5 / picosecond
MDstepsize = 0.001 * picosecond
barostat = MonteCarloBarostat(1*bar, temperature)
barostat.setFrequency(0)#disabled
system.addForce(barostat)
integrator = LangevinIntegrator(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond)

#sets up platform
platform_name = 'OpenCL'
platform = Platform.getPlatformByName(platform_name)
properties = {}

simulation = Simulation(prmtop.topology, system, integrator,platform, properties)
print ("Using platform %s" % simulation.context.getPlatform().getName())
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

state = simulation.context.getState(getEnergy = True)
print("Potential energy before minimization =", simulation.context.getState(getEnergy = True).getPotentialEnergy())

print("Running another minimization with user restraints...")
simulation.minimizeEnergy()

print("Equilibration ...")

totalSteps = 50000
steps_per_cycle = 5000
number_of_cycles = int(totalSteps/steps_per_cycle)
simulation.reporters.append(StateDataReporter(stdout, steps_per_cycle, step=True, potentialEnergy = True, temperature=True, volume=True))
simulation.reporters.append(DCDReporter(jobname + "_equil.dcd", steps_per_cycle))

simulation.step(totalSteps)

#saves checkpoint
print( "SaveState ...")
simulation.saveState(jobname + '_equil.xml')

#saves a pdb file
positions = simulation.context.getState(getPositions=True).getPositions()
boxsize = simulation.context.getState().getPeriodicBoxVectors()
simulation.topology.setPeriodicBoxVectors(boxsize)
with open(jobname + '_equil.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)
    
end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
