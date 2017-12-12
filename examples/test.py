from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os, re,time, shutil
from dmsreader import *
from datetime import datetime
from BEDAMplugin import *

lambdaList = [@lambda@]
print("Started at: " + str(time.asctime()))
start=datetime.now()
binding_file = 'test_@n@.out'
f = open(binding_file, 'w')

print("lambda = %f" %(@lambda@))
shutil.copyfile('test_rcpt_@nm1@.dms','test_rcpt_@n@.dms') 
shutil.copyfile('test_lig_@nm1@.dms','test_lig_@n@.dms') 

testDes = DesmondDMSFile('test_rcpt_@n@.dms','test_lig_@n@.dms',BEDAM=True) 

system = testDes.createSystem(nonbondedMethod=CutoffNonPeriodic,nonbondedCutoff=2.0*nanometer, OPLS = True,implicitSolvent='@implicitsolvent@')

temperature = @temperature@ * kelvin
frictionCoeff = 0.5 / picosecond
MDstepsize = 0.002 * picosecond
natoms_ligand = 12
rcpt_atom_restr = 109
lig_atom_restr = 154
kf = 1225.0
r0 = 0.5


integrator = LangevinIntegratorBEDAM(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond,natoms_ligand,@lambda@,rcpt_atom_restr,lig_atom_restr,kf, r0)

platform = Platform.getPlatformByName('@platform@')

properties = {}

if '@platform@'=='OpenCL':
    #expected "platformid:deviceid" or empty
    device = "@pn@"
    m = re.match("(\d+):(\d+)", device)
    if m:
        platformid = m.group(1)
        deviceid = m.group(2)
        properties["OpenCLPlatformIndex"] = platformid
        properties["DeviceIndex"] = deviceid
        print("Using platform id: %s, device id: %s" % ( platformid ,  deviceid) )

simulation = Simulation(testDes.topology, system, integrator,platform, properties)
simulation.context.setPositions(testDes.positions)
simulation.context.setVelocities(testDes.velocities)
testDes.get_orig_force_parameters(system)
context=simulation.context
"""
#calculate the binding energy before MD simulation
print "calculating binding energy"
testDes.binding_energy_calculation(system,context,None)
state=context.getState(getEnergy=True)
U1=state.getPotentialEnergy()
print(state.getPotentialEnergy())
testDes.set_orig_force_parameters(system,context)
testDes.binding_energy_calculation(system,context,True)
state=context.getState(getEnergy=True)
U0=state.getPotentialEnergy()
print(state.getPotentialEnergy())
print("initial binding energy="+str(U1-U0))
#ended
"""    
state = simulation.context.getState(getEnergy = True,getForces=True)
print "Using platform %s" % simulation.context.getPlatform().getName()
print(state.getPotentialEnergy())
simulation.minimizeEnergy()

stepId = @stepgap@
totalSteps = @totalsteps@
loopStep = totalSteps/stepId

simulation.reporters.append(StateDataReporter(stdout, stepId, step=True, temperature=True))
#simulation.reporters.append(DCDReporter('bcd_benzene.dcd', stepId))
#simulation.step(totalStep)

print "binding energy calculation starting"    
for i in range(loopStep):
    simulation.step(stepId)
    if simulation.currentStep%stepId==0:
        print "calculating binding energy"
        testDes.binding_energy_calculation(system,context,None)
        state=context.getState(getEnergy=True)
        U1=state.getPotentialEnergy()
	Uk=state.getKineticEnergy() #kinect energy
	Ut = (U1._value+Uk._value)/4.18 #total energy
        print(state.getPotentialEnergy())
        testDes.set_orig_force_parameters(system,context)
        testDes.binding_energy_calculation(system,context,True)
        state=context.getState(getEnergy=True)
        U0=state.getPotentialEnergy()
        print(state.getPotentialEnergy())
	Ub = (U1._value-U0._value)/4.18 #kcal/mol instead of KJ/mol
        print("binding energy="+str(U1-U0))
        #if i>10:	
        f.write("%i %f %f %f %f\n" % (simulation.currentStep, @temperature@,@lambda@, Ub, Ut))	
        testDes.set_orig_force_parameters(system,context)
f.close()
positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()
print "Updating positions and velocities"
testDes.setPositions(positions)
testDes.setVelocities(velocities)
testDes.close()

end=datetime.now()
elapsed=end - start
print("elapsed time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
