from __future__ import print_function

from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import os, re,time, shutil, math
from simtk.openmm.app.desmonddmsfile import *
from datetime import datetime

print("Started at: " + str(time.asctime()))
start=datetime.now()

output_file = 't4l_@n@.out'
f = open(output_file, 'w')


temperature = @temperature@ * kelvin

print("temperature = ", temperature)
input_structure_file  = 't4l_@nm1@.dms'
tmp_structure_file = 't4l_tmp.dms'
output_structure_file  = 't4l_@n@.dms'

shutil.copyfile(input_structure_file, tmp_structure_file)

testDes = DesmondDMSFile([tmp_structure_file]) 
# TODO: generalize cutoff treatment
system = testDes.createSystem(nonbondedMethod=NoCutoff, OPLS = True, implicitSolvent= 'AGBNP')

platform_name = '@platform@'
platform = Platform.getPlatformByName(platform_name)

properties = {}

if platform_name =='OpenCL':
    #expected "platformid:deviceid" or empty
    device = "@pn@"
    m = re.match("(\d+):(\d+)", device)
    if m:
        platformid = m.group(1)
        deviceid = m.group(2)
        properties["OpenCLPlatformIndex"] = platformid
        properties["DeviceIndex"] = deviceid
        print("Using platform id: %s, device id: %s" % ( platformid ,  deviceid) )

frictionCoeff = 0.5 / picosecond
MDstepsize = 0.001 * picosecond

integrator = LangevinIntegrator(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond)

simulation = Simulation(testDes.topology, system, integrator,platform, properties)
simulation.context.setPositions(testDes.positions)
simulation.context.setVelocities(testDes.velocities)

totalSteps = 5000
nprnt = 1000
ntrj = 1000
simulation.reporters.append(StateDataReporter(stdout, nprnt, step=True, temperature=True))
simulation.reporters.append(DCDReporter("t4l_@n@.dcd", ntrj))
simulation.reporters.append(PDBReporter("t4l_@n@.pdb", totalSteps))

loops = totalSteps/nprnt
start=datetime.now()
step = 0
for i in range(loops):
    simulation.step(nprnt)
    pot_energy = simulation.context.getState(getEnergy = True).getPotentialEnergy()
    f.write("%f %f\n" % (temperature/kelvin, pot_energy/kilocalorie_per_mole))
    f.flush()
    step += nprnt
end=datetime.now()

positions = simulation.context.getState(getPositions=True).getPositions()
velocities = simulation.context.getState(getVelocities=True).getVelocities()
testDes.setPositions(positions)
testDes.setVelocities(velocities)
testDes.close()

f.close()

shutil.copyfile(tmp_structure_file, output_structure_file)

elapsed=end - start
print("MD time="+str(elapsed.seconds+elapsed.microseconds*1e-6)+"s")
