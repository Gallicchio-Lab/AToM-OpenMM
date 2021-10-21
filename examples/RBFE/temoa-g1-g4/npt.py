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

jobname = "temoa-g1-g4"

#define the receptor atoms and the ligands atoms based on their resids
rcpt_resid = 1
lig1_resid = 2
lig2_resid = 3

#displacement vector and Vsite offsets
displ = [ 22.0, 22.0, 22.0 ]
displacement      = [  displ[i] for i in range(3) ] * angstrom
lig1_restr_offset = [  0.       for i in range(3) ] * angstrom
lig2_restr_offset = [  displ[i] for i in range(3) ] * angstrom

#reference atoms of the ligands for the alignment restraints
#these are molecular indexes (0 is the first atom of each ligand)
refatoms_lig1 = [8, 6, 4]
refatoms_lig2 = [3, 5, 1]

#defined the thermodynamic/alchemical state
#the system is prepared at the alchemical intermediate state at lambda=1/2
temperature = 300.0 * kelvin
lmbd = 0.50
lambda1 = lmbd
lambda2 = lmbd
alpha = 0.0 / kilocalorie_per_mole
u0 = 0.0 * kilocalorie_per_mole
w0coeff = 0.0 * kilocalorie_per_mole
umsc =  100.0 * kilocalorie_per_mole
ubcore = 50.0 * kilocalorie_per_mole
acore = 0.062500

#load system
prmtop = AmberPrmtopFile(jobname + '.prmtop')
inpcrd = AmberInpcrdFile(jobname + '.inpcrd')
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                             constraints=HBonds)

#load the ATM Meta Force facility. Among other things the initializer
#sorts Forces into groups 
atm_utils = ATMMetaForceUtils(system)

number_of_atoms = prmtop.topology.getNumAtoms()

#list of receptor atoms and ligands
rcpt_atoms = []
for at in prmtop.topology.atoms():
    if int(at.residue.id) == rcpt_resid:
        rcpt_atoms.append(at.index)
        
lig1_atoms = []
for at in prmtop.topology.atoms():
    if int(at.residue.id) == lig1_resid:
        lig1_atoms.append(at.index)
        
lig2_atoms = []
for at in prmtop.topology.atoms():
    if int(at.residue.id) == lig2_resid:
        lig2_atoms.append(at.index)

#Vsite restraints for each ligand. 
#the CMs are defined by all atoms of each molecule
rcpt_cm_atoms = rcpt_atoms
lig1_cm_atoms = lig1_atoms
lig2_cm_atoms = lig2_atoms
#Output lists of atoms so they can be placed into the asyncre cntl file if needed
print("List of ligand 1 CM atoms (LIGAND1_CM_ATOMS):")
print(lig1_cm_atoms)
print("List of ligand 2 CM atoms (LIGAND2_CM_ATOMS):")
print(lig2_cm_atoms)
print("List of receptor CM atoms (RCPT_CM_ATOMS):")
print(rcpt_cm_atoms)
#force constant and tolerance     
kf = 25.0 * kilocalorie_per_mole/angstrom**2 #force constant for Vsite CM-CM restraint
r0 = 5 * angstrom #radius of Vsite sphere
#Vsite restraint for the first ligand, assumed in the binding site
atm_utils.addRestraintForce(lig_cm_particles = lig1_cm_atoms,
                            rcpt_cm_particles = rcpt_cm_atoms,
                            kfcm = kf,
                            tolcm = r0,
                            offset = lig1_restr_offset)
#Vsite restraint for the second ligand, assumed in the solvent bulk
atm_utils.addRestraintForce(lig_cm_particles = lig2_cm_atoms,
                            rcpt_cm_particles = rcpt_cm_atoms,
                            kfcm = kf,
                            tolcm = r0,
                            offset = lig2_restr_offset)

#alignment restraints, to force the two ligands to remains rougly superimposed
#when they are both translated into the binding site
lig1_ref_atoms  = [ refatoms_lig1[i]+lig1_atoms[0] for i in range(3)]
lig2_ref_atoms  = [ refatoms_lig2[i]+lig2_atoms[0] for i in range(3)]
atm_utils.addAlignmentForce(liga_ref_particles = lig1_ref_atoms,
                            ligb_ref_particles = lig2_ref_atoms,
                            kfdispl = 2.5 * kilocalorie_per_mole/angstrom**2,
                            ktheta =  10.0 * kilocalorie_per_mole,
                            kpsi =  10.0 * kilocalorie_per_mole,
                            offset = lig2_restr_offset)

#restrain all heavy atoms of receptor and ligands to equilibrate only the solvent
fc = 25.0 * kilocalorie_per_mole/angstrom**2
tol = 0.5 * angstrom
hydrogen = Element.getByAtomicNumber(1)
posrestr_atoms = []
for at in prmtop.topology.atoms():
    if ( int(at.residue.id) == rcpt_resid or int(at.residue.id) == lig1_resid or int(at.residue.id) == lig2_resid ) and at.element is not hydrogen:
        posrestr_atoms.append(at.index)
atm_utils.addPosRestraints(posrestr_atoms, inpcrd.positions, fc, tol)

#create ATM Force
atmforce = ATMMetaForce(lambda1, lambda2,  alpha * kilojoules_per_mole, u0/kilojoules_per_mole, w0coeff/kilojoules_per_mole, umsc/kilojoules_per_mole, ubcore/kilojoules_per_mole, acore )
#adds all atoms to the force with zero displacement
for at in prmtop.topology.atoms():
    atmforce.addParticle(int(at.id)-1, 0., 0., 0.)
#the ligand atoms get displaced, ligand 1 from binding site to the solvent bulk
#and ligand 2 from the bulk solvent to the binding site
for i in lig1_atoms:
    atmforce.setParticleParameters(i, i, displ[0] * angstrom, displ[1] * angstrom, displ[2] * angstrom)
for i in lig2_atoms:
    atmforce.setParticleParameters(i, i, -displ[0] * angstrom, -displ[1] * angstrom, -displ[2] * angstrom)
#The ATMMetaForce assumes to be in force group 3, it looks for bonded forces in group 1 and non-bonded forces in group 2
#Bonded forces are those that are not expected to change when the ligand is displaced.
#Conversely, non-bonded forces change when the ligand is displaced. 
#The ATMMetaForceUtils() initializer (above) automatically sorts the system forces in groups 1 except for non-bonded forces
#that are placed in group 2.
atmforce.setForceGroup(3)
system.addForce(atmforce)


#Set up Langevin integrator with NPT barostat
temperature = 300 * kelvin
frictionCoeff = 0.5 / picosecond
MDstepsize = 0.001 * picosecond
barostat = MonteCarloBarostat(1*bar, temperature)
barostat.setForceGroup(1)
system.addForce(barostat)
integrator = LangevinIntegrator(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond)
#MD is conducted using forces from groups 1 and 3 only. Group 1 are bonded forces that are calculated once.
#Group 3 contains the ATMMetaForce that computes the non-bonded forces before and after the ligand is displaced and
#it then combines them according to the alchemical potential.
integrator.setIntegrationForceGroups({1,3})

platform_name = 'OpenCL'
platform = Platform.getPlatformByName(platform_name)
properties = {}

simulation = Simulation(prmtop.topology, system, integrator,platform, properties)
print ("Using platform %s" % simulation.context.getPlatform().getName())
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

pote = simulation.context.getState(getEnergy = True, groups = {1,3}).getPotentialEnergy()

print( "LoadState ...")
simulation.loadState(jobname + '_mintherm.xml')

print("Potential Energy =", simulation.context.getState(getEnergy = True, groups = {1,3}).getPotentialEnergy())

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
