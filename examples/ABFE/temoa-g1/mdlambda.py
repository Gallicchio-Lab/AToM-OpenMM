from __future__ import print_function

from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os, re,time, shutil, math
from datetime import datetime

from atmmetaforce import *

#
# This jobs performs annealing of the system from the bound state (lambda = 0) to the
# symmetric alchemical intermediate (lambda = 1/2) 
#

print("Started at: " + str(time.asctime()))
start=datetime.now()

jobname = "temoa-g1"

#temperature and initial alchemical parameters. lambda=0 is the bound state.
temperature = 300.0 * kelvin
lmbd = 0.0
lambda1 = lmbd
lambda2 = lmbd
alpha = 0.0 / kilocalorie_per_mole
u0 = 0.0 * kilocalorie_per_mole
w0coeff = 0.0 * kilocalorie_per_mole
umsc =  1000.0 * kilocalorie_per_mole
ubcore = 500.0 * kilocalorie_per_mole
acore = 0.062500

rcpt_resid = 1
lig_resid = 2

#displacement vector
displ = [ 22.0, 22.0, 22.0 ]
displacement      = [  displ[i] for i in range(3) ] * angstrom

#this is the offset of the binding site restraint relative to the host. See below.
#In leg1 the restraint is centered on the host so the offset is zero
lig_restr_offset = [  0.       for i in range(3) ] * angstrom


prmtop = AmberPrmtopFile(jobname + '.prmtop')
inpcrd = AmberInpcrdFile(jobname + '.inpcrd')
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.9*nanometer,
                             constraints=HBonds)
atm_utils = ATMMetaForceUtils(system)

number_of_atoms = prmtop.topology.getNumAtoms()

rcpt_atoms = []
for at in prmtop.topology.atoms():
    if int(at.residue.id) == rcpt_resid:
        rcpt_atoms.append(at.index)
        
lig_atoms = []
for at in prmtop.topology.atoms():
    if int(at.residue.id) == lig_resid:
        lig_atoms.append(at.index)
print("List of ligand atoms (LIGAND_ATOMS):")
print(lig_atoms)

#these are the atoms of the receptor and the ligand that define their centroids
lig_cm_atoms = lig_atoms
print("List of ligand CM atoms (LIGAND_CM_ATOMS):")
print(lig_cm_atoms)

rcpt_cm_atoms = rcpt_atoms
print("List of receptor CM atoms (RCPT_CM_ATOMS):")
print(rcpt_cm_atoms)


#define the binding site restraint. It is a flat-bottom central potential based on the CM-CM distance
#between the receptor and the ligand. r0 is the tolerance and kf is the force constant
#of the quadratic penalty term if the tolerance is exceeded. r0 is also regarded the radius of the
#of the spherical binding site region of volume Vsite = (4/3) pi r0^3
kf = 25.0 * kilocalorie_per_mole/angstrom**2
r0 = 5 * angstrom
atm_utils.addRestraintForce(lig_cm_particles = lig_cm_atoms,
                            rcpt_cm_particles = rcpt_cm_atoms,
                            kfcm = kf,
                            tolcm = r0,
                            offset = lig_restr_offset)

#receptor positional restraints
fc = 25.0 * kilocalorie_per_mole/angstrom**2
tol = 0.5 * angstrom
carbon = re.compile("^C.*")
posrestr_atoms = []
for at in prmtop.topology.atoms():
    if int(at.residue.id) == rcpt_resid and carbon.match(at.name) and at.index < 40:
        posrestr_atoms.append(at.index)
atm_utils.addPosRestraints(posrestr_atoms, inpcrd.positions, fc, tol)
print("List of atoms whose positions are restrained (POS_RESTRAINED_ATOMS):")
print(posrestr_atoms)
        
#create ATM Force
atmforce = ATMMetaForce(lambda1, lambda2,  alpha * kilojoules_per_mole, u0/kilojoules_per_mole, w0coeff/kilojoules_per_mole, umsc/kilojoules_per_mole, ubcore/kilojoules_per_mole, acore )
#adds atoms to the force with zero displacement
for at in prmtop.topology.atoms():
    atmforce.addParticle(at.index, 0., 0., 0.)
#the ligand atoms get displaced
for i in lig_atoms:
    atmforce.setParticleParameters(i, i, displ[0] * angstrom, displ[1] * angstrom, displ[2] * angstrom)
#The ATMMetaForce assumes to be in force group 3, it looks for bonded forces in group 1 and non-bonded forces in group 2
#Bonded forces are those that are not expected to change when the ligand is displaced.
#Conversely, non-bonded forces change when the ligand is displaced. 
#The ATMMetaForceUtils() initializer (above) automatically sorts the system forces in groups 1 except for non-bonded forces
#that are placed in group 2.
atmforce.setForceGroup(3)
system.addForce(atmforce)

#add barostat, but disables it to run NVT.
#it is assigned to force group 1
barostat = MonteCarloBarostat(1*bar, temperature)
barostat.setForceGroup(1)
barostat.setFrequency(0)#disabled
system.addForce(barostat)

#Set up Langevin integrator
frictionCoeff = 0.5 / picosecond
MDstepsize = 0.001 * picosecond
#MD is conducted using forces from groups 1 and 3 only. Group 1 are bonded forces that are calculated once.
#Group 3 contains the ATMMetaForce that computes the non-bonded forces before and after the ligand is displaced and
#it then combines them according to the alchemical potential.
integrator = MTSLangevinIntegrator(temperature, frictionCoeff, MDstepsize, [(3,1), (1,1)])
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

#For reasons unknown an initial energy calculation appears to be required before MD
pote = simulation.context.getState(getEnergy = True, groups = {1,3}).getPotentialEnergy()

print( "LoadState ...")
simulation.loadState(jobname + '_equil.xml')

#sets the ATM alchemical parameters
simulation.context.setParameter(atmforce.Lambda1(), lambda1)
simulation.context.setParameter(atmforce.Lambda2(), lambda2)
simulation.context.setParameter(atmforce.Alpha(), alpha *kilojoules_per_mole)
simulation.context.setParameter(atmforce.U0(), u0 /kilojoules_per_mole)
simulation.context.setParameter(atmforce.W0(), w0coeff /kilojoules_per_mole)
simulation.context.setParameter(atmforce.Umax(), umsc /kilojoules_per_mole)
simulation.context.setParameter(atmforce.Ubcore(), ubcore /kilojoules_per_mole)
simulation.context.setParameter(atmforce.Acore(), acore)

print("Potential Energy =", simulation.context.getState(getEnergy = True, groups = {1,3}).getPotentialEnergy())

print("Annealing to lambda = 1/2 ...")

#goes from lambda=0 to lambda=1/2 in 250,000 MD steps updating lambda every 1000 steps
totalSteps = 250000
steps_per_cycle = 1000
number_of_cycles = int(totalSteps/steps_per_cycle)
deltalambda = (0.5 - 0.0)/float(number_of_cycles)
simulation.reporters.append(StateDataReporter(stdout, steps_per_cycle, step=True, potentialEnergy = True, temperature=True))
simulation.reporters.append(DCDReporter(jobname + "_mdlambda.dcd", steps_per_cycle))

#binding energy values and other parameters are recorded in this file
f = open(jobname + "_mdlambda.out", 'w')

for i in range(number_of_cycles):
    simulation.step(steps_per_cycle)
    state = simulation.context.getState(getEnergy = True, groups = {1,3})
    pot_energy = (state.getPotentialEnergy()).value_in_unit(kilocalorie_per_mole)
    pert_energy = (atmforce.getPerturbationEnergy(simulation.context)).value_in_unit(kilocalorie_per_mole)
    l1 = simulation.context.getParameter(atmforce.Lambda1())
    l2 = simulation.context.getParameter(atmforce.Lambda2())
    a = simulation.context.getParameter(atmforce.Alpha()) / kilojoules_per_mole
    umid = simulation.context.getParameter(atmforce.U0()) * kilojoules_per_mole
    w0 = simulation.context.getParameter(atmforce.W0()) * kilojoules_per_mole
    print("%f %f %f %f %f %f %f %f %f" % (temperature/kelvin,lmbd, l1, l2, a*kilocalorie_per_mole, umid/kilocalorie_per_mole, w0/kilocalorie_per_mole, pot_energy, pert_energy), file=f )
    f.flush()
    lmbd += deltalambda
    #for a linear alchemical potential lambda1 = lambda2 = lambda
    lambda1 += deltalambda
    lambda2 += deltalambda
    simulation.context.setParameter(atmforce.Lambda1(), lambda1)
    simulation.context.setParameter(atmforce.Lambda2(), lambda2)

print( "SaveState ...")

#saves a checkpoint file and a pdb file at lambda=1/2.
#for historical reasons it has a "0" in the file name
simulation.saveState(jobname + '_0.xml')
positions = simulation.context.getState(getPositions=True).getPositions()
boxsize = simulation.context.getState().getPeriodicBoxVectors()
simulation.topology.setPeriodicBoxVectors(boxsize)
with open(jobname + '_0.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)
  
#create a checkpoint file and a pdb file the ligand displaced in the bulk for leg2
positions = simulation.context.getState(getPositions=True).getPositions()
for i in lig_atoms:
    positions[i] += displacement
simulation.context.setPositions(positions)
simulation.saveState(jobname + '_0_displaced.xml')
with open(jobname + '_0_displaced.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)

end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
