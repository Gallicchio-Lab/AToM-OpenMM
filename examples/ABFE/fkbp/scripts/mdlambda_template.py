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

displ = [ <DISPLX>, <DISPLY>, <DISPLZ> ]
off = [ <OFFX>, <OFFY>, <OFFZ> ]
displacement      = [  displ[i] for i in range(3) ] * angstrom
lig_restr_offset = [  off[i]       for i in range(3) ] * angstrom


lig_atoms = [<LIGATOMS>]
rcpt_cm_atoms = [<VSITERECEPTORATOMS>]
restrained_atoms = [<RESTRAINEDATOMS>]

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

prmtop = AmberPrmtopFile(jobname + '.prmtop')
inpcrd = AmberInpcrdFile(jobname + '.inpcrd')
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                             constraints=HBonds)
atm_utils = ATMMetaForceUtils(system)


#rcpt_atoms = []
#rcpt_atoms = [at.index for at in prmtop.topology.atoms() if at.residue.name =='PA' or at.residue.name == 'PC' or at.residue.name == 'OL' or at.residue.name == 'CHL' ]

"""
kf = 25.0 * kilocalorie_per_mole/angstrom**2 #force constant for Vsite CM-CM restraint
r0 = 5 * angstrom #radius of Vsite sphere
offz =  lig_restr_offset[2]
zrestrforce = mm.CustomCentroidBondForce(2,"0.5*kf*( step(d)*max(0,d-r0)^2 + step(-d)*max(0,-d-r0)^2 ) ; d = z2 - offz - z1")
system.addForce(zrestrforce)
zrestrforce.addPerBondParameter("kf")
zrestrforce.addPerBondParameter("r0")
zrestrforce.addPerBondParameter("offz")
zrestrforce.setForceGroup(1)
zrestrforce.addGroup(rcpt_atoms)
zrestrforce.addGroup(lig_atoms)
zrestrforce.addBond([0,1], [kf, r0, offz ])
"""


lig_cm_atoms = lig_atoms

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

#create ATM Force
atmforce = ATMMetaForce(lambda1, lambda2,  alpha * kilojoules_per_mole, u0/kilojoules_per_mole, w0coeff/kilojoules_per_mole, umsc/kilojoules_per_mole, ubcore/kilojoules_per_mole, acore )

for at in prmtop.topology.atoms():
    atmforce.addParticle(at.index, 0., 0., 0.)

for i in lig_atoms:
    atmforce.setParticleParameters(i, i, displ[0] * angstrom, displ[1] * angstrom, displ[2] * angstrom)

atmforce.setForceGroup(3)

system.addForce(atmforce)

#add barostat
barostat = MonteCarloBarostat(1*bar, temperature)
barostat.setForceGroup(1)
barostat.setFrequency(0)#disabled
system.addForce(barostat)

temperature = 300 * kelvin
frictionCoeff = 0.5 / picosecond
MDstepsize = 0.001 * picosecond

integrator = LangevinIntegrator(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond)
integrator.setIntegrationForceGroups({1,3})

platform_name = 'OpenCL'
#platform_name = 'Reference'
platform = Platform.getPlatformByName(platform_name)

properties = {}

simulation = Simulation(prmtop.topology, system, integrator,platform, properties)
print ("Using platform %s" % simulation.context.getPlatform().getName())
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

state = simulation.context.getState(getEnergy = True, groups = {1,3})
pote = state.getPotentialEnergy()

print( "LoadState ...")
simulation.loadState(jobname + '_equil.xml')

#override ATM parameters from checkpoint file
simulation.context.setParameter(atmforce.Lambda1(), lambda1)
simulation.context.setParameter(atmforce.Lambda2(), lambda2)
simulation.context.setParameter(atmforce.Alpha(), alpha *kilojoules_per_mole)
simulation.context.setParameter(atmforce.U0(), u0 /kilojoules_per_mole)
simulation.context.setParameter(atmforce.W0(), w0coeff /kilojoules_per_mole)
simulation.context.setParameter(atmforce.Umax(), umsc /kilojoules_per_mole)
simulation.context.setParameter(atmforce.Ubcore(), ubcore /kilojoules_per_mole)
simulation.context.setParameter(atmforce.Acore(), acore)

state = simulation.context.getState(getEnergy = True, groups = {1,3})
print("Potential Energy =", state.getPotentialEnergy())

print("Annealing to lambda = 1/2 ...")

stepId = 1000
totalSteps = 250000
loopStep = int(totalSteps/stepId)
simulation.reporters.append(StateDataReporter(stdout, stepId, step=True, potentialEnergy = True, temperature=True))
simulation.reporters.append(DCDReporter(jobname + ".dcd", stepId))

binding_file = jobname + '.out'
f = open(binding_file, 'w')

deltalambda = (0.5 - 0.0)/float(loopStep)

for i in range(loopStep):
    simulation.step(stepId)
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
    lambda1 += deltalambda
    lambda2 += deltalambda
    simulation.context.setParameter(atmforce.Lambda1(), lambda1)
    simulation.context.setParameter(atmforce.Lambda2(), lambda2)

print( "SaveState ...")
simulation.saveState(jobname + "_0.xml")

#create a checkpoint with the ligands swapped for leg2
positions = simulation.context.getState(getPositions=True).getPositions()
for i in lig_atoms:
    positions[i] += displacement
simulation.context.setPositions(positions)
simulation.saveState(jobname + '_0_displaced.xml')

#save a pdb file with the displaced ligand for visualization
positions = simulation.context.getState(getPositions=True).getPositions()
with open(jobname + '_0_displaced.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)

end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
