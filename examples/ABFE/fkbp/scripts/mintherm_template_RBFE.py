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
displacement      = [  displ[i] for i in range(3) ] * angstrom
lig1_restr_offset = [  0.       for i in range(3) ] * angstrom
lig2_restr_offset = [  displ[i] for i in range(3) ] * angstrom

lig1_atoms = [ <LIG1ATOMS> ]
lig2_atoms = [ <LIG2ATOMS> ]
refatoms_lig1 = [ <REFERENCEATOMS1> ]
refatoms_lig2 = [ <REFERENCEATOMS2> ]
rcpt_cm_atoms = [ <VSITERECEPTORATOMS> ]

#define the thermodynamic/alchemical state
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

#Vsite restraints
lig1_cm_atoms = lig1_atoms
lig2_cm_atoms = lig2_atoms
kf = 25.0 * kilocalorie_per_mole/angstrom**2 #force constant for Vsite CM-CM restraint
r0 = 5 * angstrom #radius of Vsite sphere
atm_utils.addRestraintForce(lig_cm_particles = lig1_cm_atoms,
                            rcpt_cm_particles = rcpt_cm_atoms,
                            kfcm = kf,
                            tolcm = r0,
                            offset = lig1_restr_offset)

atm_utils.addRestraintForce(lig_cm_particles = lig2_cm_atoms,
                            rcpt_cm_particles = rcpt_cm_atoms,
                            kfcm = kf,
                            tolcm = r0,
                            offset = lig2_restr_offset)

#alignment restraint
kfdispl = 2.5 * kilocalorie_per_mole/angstrom**2
ktheta = 10.0 * kilocalorie_per_mole
kpsi   = 10.0 * kilocalorie_per_mole
lig1_ref_atoms  = [ refatoms_lig1[i]+lig1_atoms[0] for i in range(3)]
lig2_ref_atoms  = [ refatoms_lig2[i]+lig2_atoms[0] for i in range(3)]
atm_utils.addAlignmentForce(liga_ref_particles = lig1_ref_atoms,
                            ligb_ref_particles = lig2_ref_atoms,
                            kfdispl = 2.5 * kilocalorie_per_mole/angstrom**2,
                            ktheta =  10.0 * kilocalorie_per_mole,
                            kpsi =  10.0 * kilocalorie_per_mole,
                            offset = lig2_restr_offset)

#restrain all heavy atoms of receptor and ligands to equilibrate only the solvent
#this assumes that ligand 2 is the last residue of the receptor-ligand system
fc = 25.0 * kilocalorie_per_mole/angstrom**2
tol = 0.5 * angstrom
hydrogen = Element.getByAtomicNumber(1)
posrestr_atoms = []
last_not_solvent_atom = lig2_atoms[len(lig2_atoms)-1]
for at in prmtop.topology.atoms():
    if at.index <= last_not_solvent_atom:
        posrestr_atoms.append(at.index)
atm_utils.addPosRestraints(posrestr_atoms, inpcrd.positions, fc, tol)

#create ATM Force
atmforce = ATMMetaForce(lambda1, lambda2,  alpha * kilojoules_per_mole, u0/kilojoules_per_mole, w0coeff/kilojoules_per_mole, umsc/kilojoules_per_mole, ubcore/kilojoules_per_mole, acore )
#adds all atoms to the force with zero displacement
for at in prmtop.topology.atoms():
    atmforce.addParticle(at.index, 0., 0., 0.)
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

#Set up Langevin integrator
initial_temperature = 50 * kelvin
final_temperature = 300 * kelvin
temperature = initial_temperature
frictionCoeff = 0.5 / picosecond
MDstepsize = 0.001 * picosecond
barostat = MonteCarloBarostat(1*bar, temperature)
barostat.setForceGroup(1)
saved_barostat_frequency = barostat.getFrequency()
barostat.setFrequency(0)#disabled
system.addForce(barostat)
integrator = LangevinIntegrator(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond)
#MD is conducted using forces from groups 1 and 3 only. Group 1 are bonded forces that are calculated once.
#Group 3 contains the ATMMetaForce that computes the non-bonded forces before and after the ligand is displaced and
#it then combines them according to the alchemical potential.
integrator.setIntegrationForceGroups({1,3})

#sets up platform
platform_name = 'OpenCL'
platform = Platform.getPlatformByName(platform_name)
properties = {}

simulation = Simulation(prmtop.topology, system, integrator,platform, properties)
print ("Using platform %s" % simulation.context.getPlatform().getName())
simulation.context.setPositions(inpcrd.positions)
if inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

state = simulation.context.getState(getEnergy = True, groups = {1,3})
print("Potential energy before minimization =", simulation.context.getState(getEnergy = True, groups = {1,3}).getPotentialEnergy())

print("Energy minimizing the system ...")
simulation.minimizeEnergy()

print("Potential energy after minimization =", simulation.context.getState(getEnergy = True, groups = {1,3}).getPotentialEnergy())
  
print("Thermalization ...")

totalSteps = 50000
steps_per_cycle = 5000
number_of_cycles = int(totalSteps/steps_per_cycle)
delta_temperature = (final_temperature - initial_temperature)/number_of_cycles
simulation.reporters.append(StateDataReporter(stdout, steps_per_cycle, step=True, potentialEnergy = True, temperature=True, volume=True))

#binding energy values and other parameters are recorded in this file
f = open(jobname + "_mintherm.out", 'w')

#MD with temperature ramp
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
    #prepare system for new temperature
    temperature = temperature + delta_temperature
    integrator.setTemperature(temperature)
    
print("NPT equilibration ...")
barostat.setFrequency(saved_barostat_frequency)#enabled

#MD at constant pressure
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

#saves checkpoint
print( "SaveState ...")
simulation.saveState(jobname + '_mintherm.xml')

#saves a pdb file
positions = simulation.context.getState(getPositions=True).getPositions()
boxsize = simulation.context.getState().getPeriodicBoxVectors()
simulation.topology.setPeriodicBoxVectors(boxsize)
with open(jobname + '_mintherm.pdb', 'w') as output:
  PDBFile.writeFile(simulation.topology, positions, output)
    
end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
