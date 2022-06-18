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

temperature = 300.0 * kelvin

rcpt_resid = 1
lig_resid = 2

prmtop = AmberPrmtopFile(jobname + '.prmtop')
inpcrd = AmberInpcrdFile(jobname + '.inpcrd')
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.9*nanometer,
                             constraints=HBonds)
atm_utils = ATMMetaForceUtils(system)

#restrain the positions all heavy atoms of receptor and ligand to relax only the solvent
fc = 25.0 * kilocalorie_per_mole/angstrom**2
tol = 0.5 * angstrom
hydrogen = Element.getByAtomicNumber(1)
posrestr_atoms = []
for at in prmtop.topology.atoms():
    if ( int(at.residue.id) == rcpt_resid or int(at.residue.id) == lig_resid ) and at.element is not hydrogen:
        posrestr_atoms.append(at.index)
atm_utils.addPosRestraints(posrestr_atoms, inpcrd.positions, fc, tol)

temperature = 300 * kelvin

#add barostat
barostat = MonteCarloBarostat(1*bar, temperature)
system.addForce(barostat)

#set up integrator
frictionCoeff = 0.5 / picosecond
MDstepsize = 0.001 * picosecond
nonbonded_force_group = 1
atm_utils.setNonbondedForceGroup(nonbonded_force_group)
integrator = MTSLangevinIntegrator(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond, [(0,1), (nonbonded_force_group,1)])
integrator.setConstraintTolerance(0.00001)

platform_name = 'OpenCL'
#platform_name = 'CUDA'
platform = Platform.getPlatformByName(platform_name)
properties = {}
properties["Precision"] = "mixed"

simulation = Simulation(prmtop.topology, system, integrator,platform, properties)
print ("Using platform %s" % simulation.context.getPlatform().getName())
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

pote = simulation.context.getState(getEnergy = True).getPotentialEnergy()

print( "LoadState ...")
simulation.loadState(jobname + '_mintherm.xml')

print("Potential Energy =", simulation.context.getState(getEnergy = True).getPotentialEnergy())

print("NPT ...")

nprint = 5000
totalSteps = 50000
simulation.reporters.append(StateDataReporter(stdout, nprint, step=True, potentialEnergy = True, temperature=True, volume=True))

simulation.step(totalSteps)

print( "SaveState ...")
simulation.saveState(jobname + '_npt.xml')

#save a pdb file that can be used as a topology to load .dcd files in vmd
boxsize = simulation.context.getState().getPeriodicBoxVectors()
simulation.topology.setPeriodicBoxVectors(boxsize)
positions = simulation.context.getState(getPositions=True).getPositions()
with open(jobname + '_npt.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)

end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
