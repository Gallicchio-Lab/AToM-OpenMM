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

jobname = "temoa-g1"

rcpt_resid = 1
lig_resid = 2

prmtop = AmberPrmtopFile(jobname + '.prmtop')
inpcrd = AmberInpcrdFile(jobname + '.inpcrd')
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                             constraints=HBonds)
atm_utils = ATMMetaForceUtils(system)

number_of_atoms = prmtop.topology.getNumAtoms()

#restrain the positions all heavy atoms of receptor and ligand to relax only the solvent
fc = 25.0 * kilocalorie_per_mole/angstrom**2
tol = 0.5 * angstrom
hydrogen = Element.getByAtomicNumber(1)
posrestr_atoms = []
for at in prmtop.topology.atoms():
    if ( int(at.residue.id) == rcpt_resid or int(at.residue.id) == lig_resid ) and at.element is not hydrogen:
        posrestr_atoms.append(at.index)
atm_utils.addPosRestraints(posrestr_atoms, inpcrd.positions, fc, tol)

#warm up the system from 50 K to 300 K
initial_temperature = 50 * kelvin
final_temperature = 300 * kelvin

temperature = initial_temperature
frictionCoeff = 0.5 / picosecond
MDstepsize = 0.001 * picosecond
    
integrator = LangevinIntegrator(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond)

platform_name = 'OpenCL'
platform = Platform.getPlatformByName(platform_name)

properties = {}

simulation = Simulation(prmtop.topology, system, integrator,platform, properties)

print ("Using platform %s" % simulation.context.getPlatform().getName())
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

print("Potential energy before minimization =", simulation.context.getState(getEnergy = True).getPotentialEnergy())

print("Energy minimizing the system ...")
simulation.minimizeEnergy()

print("Potential energy after minimization =", simulation.context.getState(getEnergy = True).getPotentialEnergy())
  
print("Thermalization ...")

totalSteps = 50000
steps_per_cycle = 5000
number_of_cycles = int(totalSteps/steps_per_cycle)
delta_temperature = (final_temperature - initial_temperature)/number_of_cycles
simulation.reporters.append(StateDataReporter(stdout, steps_per_cycle, step=True, potentialEnergy = True, temperature=True))

#MD with temperature ramp
for i in range(number_of_cycles):
    simulation.step(steps_per_cycle)
    temperature = temperature + delta_temperature
    integrator.setTemperature(temperature)

print( "Save state ...")
simulation.saveState(jobname + '_mintherm.xml')

#save a pdb file that can be used as a topology to load .dcd files in vmd
positions = simulation.context.getState(getPositions=True).getPositions()
boxsize = simulation.context.getState().getPeriodicBoxVectors()
simulation.topology.setPeriodicBoxVectors(boxsize)
with open(jobname + '_mintherm.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)

end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
