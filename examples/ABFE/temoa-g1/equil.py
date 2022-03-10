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

#restrain the positions of the carbon atoms of the lower cup of the host
fc = 25.0 * kilocalorie_per_mole/angstrom**2
tol = 0.5 * angstrom
carbon = re.compile("^C.*")
posrestr_atoms = []
for at in prmtop.topology.atoms():
    if int(at.residue.id) == rcpt_resid and carbon.match(at.name) and at.index < 40:
        posrestr_atoms.append(at.index)
atm_utils.addPosRestraints(posrestr_atoms, inpcrd.positions, fc, tol)

temperature = 300 * kelvin

#add barostat because the checkpoint file requires it, but disable it to run NVT
barostat = MonteCarloBarostat(1*bar, temperature)
barostat.setFrequency(0)#disabled
system.addForce(barostat)

#set up integrator
frictionCoeff = 0.5 / picosecond
MDstepsize = 0.001 * picosecond
integrator = MTSLangevinIntegrator(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond, [(1,1), (2,1)])
integrator.setConstraintTolerance(0.00001)

#platform_name = 'OpenCL'
platform_name = 'CUDA'
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
simulation.loadState(jobname + '_npt.xml')

print("Potential Energy =", simulation.context.getState(getEnergy = True).getPotentialEnergy())

print("Equilibration ...")

nprint = 5000
totalSteps = 50000
simulation.reporters.append(StateDataReporter(stdout, nprint, step=True, potentialEnergy = True, temperature=True, volume=True))

simulation.step(totalSteps)

print( "SaveState ...")
simulation.saveState(jobname + '_equil.xml')

#save a pdb file that can be used as a topology to load .dcd files in vmd
positions = simulation.context.getState(getPositions=True).getPositions()
boxsize = simulation.context.getState().getPeriodicBoxVectors()
simulation.topology.setPeriodicBoxVectors(boxsize)
with open(jobname + '_equil.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)

end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
