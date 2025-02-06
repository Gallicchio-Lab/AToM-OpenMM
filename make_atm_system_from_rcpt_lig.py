# Usage: python make_atm_rbfe_system_frompdb.py <options>
# Emilio Gallicchio, 5/2023 adapted from code by Bill Swope, 11/2021
# modified by Solmaz Azimi 6/2024

############################################
#                                          #
#   IMPORTS                                #
#                                          #
############################################

import os, sys, re
import numpy as np
from datetime import datetime
from time import time

# following for argument passing tools
import argparse

# OpenMM components
from openmm import XmlSerializer
from openmm import Vec3
from openmm.app import PDBReporter, StateDataReporter, PDBFile
from openmm.app import ForceField, Modeller
from openmm.app import PME, HBonds, NoCutoff, OBC2

# OpenFF components from the toolkit
from openff.toolkit.topology import Molecule

from openmm.unit import Quantity
from openmm.unit import angstrom, nanometer, nanometers, picoseconds, amu
from openmm import unit

from sys import stdout

# System calls - to invoke sdfTagTool from the python code
import os, sys

# OpenFF and OpenMM components for ligand force field parameters
from openff.toolkit.topology import Molecule

#from pdbfixer import PDBFixer

# OpenFF components for SystemGenerator (for input entire ssytem in pdb file)
from openmm import app
from openmmforcefields.generators import SystemGenerator

def boundingBoxSizes(positions):
    xmin = positions[0][0]
    xmax = positions[0][0]
    ymin = positions[0][1]
    ymax = positions[0][1]
    zmin = positions[0][2]
    zmax = positions[0][2]
    for i in range(len(positions)):
        x = positions[i][0]
        y = positions[i][1]
        z = positions[i][2]
        # print('type of x ', type(x))
        # print('Site ', i, ' Coord ', positions[i])
        if(x > xmax):
            xmax = x
        if(x < xmin):
            xmin = x
        if(y > ymax):
            ymax = y
        if(y < ymin):
            ymin = y
        if(z > zmax):
            zmax = z
        if(z < zmin):
            zmin = z
    return [ (xmin,xmax), (ymin,ymax), (zmin,zmax) ] 

############################################
#                                          #
#   ACTUAL CODE                            #
#                                          #
############################################

print('Generate ATM OpenMM System')
today = datetime.today()
print('\nDate and time at start: ', today.strftime('%c'))


whatItDoes = """
Produces an .xml file with the OpenMM's system for ATM binding
free energy calcuations: absolute, relative, or swapping.

The user supplies a pdb file that contains protein, cofactors, ligands, 
water, ions, probably prepared by a different software package 
(Maestro, OpenEye).  Cofactors and ligands have to be described in an
sdf file.  Then force fields are assigned to all components.

Emilio Gallicchio 5/2023
adapted from simProteinLigandWater.py script by Bill Swope 11/2021
modified by Solmaz Azimi 6/2024
"""

#####################################################
#                                                   #
#   PASS COMMAND LINE ARGUMENTS TO LOCAL VARIABLES  #
#                                                   #
#####################################################


program_start_timer = time()
parser = argparse.ArgumentParser(description=whatItDoes)

# Required input
parser.add_argument('--receptorinFile', required=True,  type=str, default=None,
                    help='Receptor in SDF (.sdf) or PDB format (.pdb)')
parser.add_argument('--LIG1SDFinFile',  required=True,  type=str, default=None,
                    help='SDF file of first ligand')
parser.add_argument('--displacement',  required=True,  type=str, default=None,
                    help='string with displacement vector in angstroms like "22. 0.0 0.0" ')
parser.add_argument('--systemXMLoutFile',  required=True,  type=str, default=None,
                    help='Name of the XML file where to save the System')
parser.add_argument('--systemPDBoutFile', required=True, type=str, default=None,
                    help='Name of the PDB file where to output to system')

# Optional input
parser.add_argument('--receptor2inFile',  required=False,  type=str, default=None,
                    help='PDB file of second receptor')
parser.add_argument('--LIG2SDFinFile',  required=False,  type=str, default=None,
                    help='SDF file of second ligand')
parser.add_argument('--cofactorsSDFFile',  required=False,  type=str, default=None,
                    help='SDF file with receptor cofactors')

parser.add_argument('--proteinForceField', required=False, type=str,
                    default='amber14-all.xml',
                    help='Force field for protein` amber14-all.xml ')
parser.add_argument('--solventForceField', required=False, type=str,
                    default='amber14/tip3p.xml',
                    help='Force field for solvent/ions` amber14/tip3p.xml ')
parser.add_argument('--ligandForceField', required=False, type=str,
                    default='openff-2.0.0',
                    help='Force field for ligand:  openff-2.0.0, gaff, or espaloma-0.3.2')
parser.add_argument('--implicitSolvent', required=False, type=str,
                    default='None',
                    help='Implicit solvent to use: HCT OBC2 GBn2. None for vacuum.')
parser.add_argument('--forcefieldJSONCachefile', required=False, type=str,
                    default=None,
                    help='Force field ligand cache database')
parser.add_argument('--hmass', required=False, type=float,
                    default=1.0,
                    help='Hydrogen mass, set it to 1.5 amu to use a 4 fs time-step')

# Arguments that are flags
parser.add_argument('--verbose', required=False, action='store_true',
                    help='Get more output with this flag')

args = vars(parser.parse_args())

#catch abfe or rbfe or lsfe
rbfe = False
if args['LIG2SDFinFile']  is not None:
    rbfe = True

lsfe = False
if args['receptor2inFile'] is not None:
    lsfe = True

# Pull data from command line into local variables
receptorfile = args['receptorinFile']
lig1sdffile = args['LIG1SDFinFile']
xmloutfile = args['systemXMLoutFile']
pdboutfile = args['systemPDBoutFile']

displ = [float(r) for r in args['displacement'].split()]
displacement = Vec3(displ[0], displ[1], displ[2]) * angstrom

if rbfe:
    lig2sdffile  = args['LIG2SDFinFile']

if lsfe:
    receptor2file = args['receptor2inFile']
    
#cofactors
cofsdffile = args['cofactorsSDFFile']
    
# force fields
proteinforcefield = args['proteinForceField']
solventforcefield = args['solventForceField']
ligandforcefield = args['ligandForceField']
ffcachefile = args['forcefieldJSONCachefile']

#implicit solvent
implsolv = args['implicitSolvent']
if implsolv == 'None':
    implsolv = None

hmass = float(args['hmass'])

# flag for printing (verbose) 
flagverbose = args['verbose']

#####################################################
#   Echo out the suppliable input parameters        #
#####################################################


print('\nUser-supplied input parameters')
print('Receptor file name:                 ', receptorfile)
if lsfe:
    print('Receptor2 file name:            ', receptor2file)
print('Protein force field:                ', proteinforcefield)
print('Solvent/ion force field             ', solventforcefield )
print('Ligand force field:                 ', ligandforcefield)
print('Ligand 1 file name:                 ', lig1sdffile)
if rbfe:
    print('Ligand 2 file name:                 ', lig2sdffile)
    print('Displacement                        ', displacement)
print('Topology PDB output file:           ', pdboutfile)
print('System XML output file:             ', xmloutfile)
print('Force field cache file:             ', ffcachefile)


print('Call ForceField for protein and water')
forcefield = ForceField(proteinforcefield,solventforcefield)
if implsolv is not None:
    if implsolv == "OBC2":
        forcefield.loadFile('implicit/obc2.xml')
    elif implsolv == "GBN2":
        forcefield.loadFile('implicit/gbn2.xml')
    elif implsolv == "HCT":
        forcefield.loadFile('implicit/hct.xml')
    elif implsolv == "Vacuum" or implsolv == "vacuum":
        pass
    else:
        print('Unknown implicit solvent %s' % implsolv)
        sys.exit(1)

# to store OpenFF molecule objects of non-protein units
ligandmolecules = []

############################################
#                                          #
#   READ AND CHARACTERIZE RECEPTORS        #
#                                          #
############################################

rcptpext = os.path.splitext(receptorfile)[1]
rcpt_ommtopology = None
rcpt_positions = None
if rcptpext == '.pdb':
    print('Receptor in PDB format')
    pdbrcpt = PDBFile(receptorfile)
    rcpt_positions = pdbrcpt.positions
    rcpt_ommtopology = pdbrcpt.topology
elif rcptpext == '.sdf':
    print('Receptor in SDF format')
    molrcpt = Molecule.from_file(receptorfile, file_format='SDF',
                                 allow_undefined_stereo=True)
    ligandmolecules.append(molrcpt)

    pos = molrcpt.conformers[0].to('angstrom').magnitude
    rcpt_positions = [Vec3(pos[i][0], pos[i][1], pos[i][2]) for i in range(pos.shape[0])] * angstrom
    
    offtopology = molrcpt.to_topology()
    rcpt_ommtopology = offtopology.to_openmm(ensure_unique_atom_names=True)
else:
    print("Error: Unrecognized receptor file name: %s" % receptorfile)

nrcpt = rcpt_ommtopology.getNumAtoms()
print('Number of atoms in receptor:', nrcpt)

print('Call Modeller: include receptor first')
modeller = Modeller(rcpt_ommtopology, rcpt_positions)

if lsfe:
    rcpt2pext = os.path.splitext(receptor2file)[1]
    rcpt2_ommtopology = None
    rcpt2_positions = None
    if rcpt2pext == '.pdb':
        print('Read second receptor:')
        pdbrcpt2 = PDBFile(receptor2file)
        rcpt2_positions = pdbrcpt2.positions
        rcpt2_ommtopology = pdbrcpt2.topology
    else:
        print("Error: Unrecognized receptor file name: %s" % receptor2file)

    nrcpt2 = rcpt2_ommtopology.getNumAtoms()
    print('Number of atoms in receptor2:', nrcpt2)

    # displace receptor
    for i in range(nrcpt2):
        rcpt2_positions[i] += displacement

    print('Call Modeller: include second receptor')
    modeller.add(rcpt2_ommtopology, rcpt2_positions)


print("Calculating receptor bounding box:")
bbox = boundingBoxSizes(rcpt_positions)
bboxsizes = [ bbox[i][1]-bbox[i][0] for i in range(3) ]
bboxfaces = [ bboxsizes[2]*bboxsizes[1], bboxsizes[2]*bboxsizes[0],  bboxsizes[1]*bboxsizes[0] ]
print("Areas of faces", bboxfaces)
smallest_direction = 0
smallest_area = bboxfaces[0]
for i in range(3):
    if bboxfaces[i] < smallest_area:
        smallest_direction = i
print("Direction of smallest area dimension:", smallest_direction)

    
############################################
#                                          #
#   READ AND CHARACTERIZE LIGANDS          #
#                                          #
############################################

if cofsdffile is not None:
    print('Read cofactor:')
    molcof = Molecule.from_file(cofsdffile, file_format='SDF',
                             allow_undefined_stereo=True)
    ligandmolecules.append(molcof)
    molcof_ommtopology = molcof.to_topology().to_openmm(ensure_unique_atom_names=True)

    #assign the residue name, assumes one residue
    resfile = os.path.split(cofsdffile)[1]
    resname = os.path.splitext(resfile)[0]
    for residue in molcof_ommtopology.residues():
        residue.name = resname.upper()

    pos = molcof.conformers[0].to('angstrom').magnitude
    molcof_positions = [Vec3(pos[i][0], pos[i][1], pos[i][2]) for i in range(pos.shape[0])] * angstrom
    ncof = molcof_ommtopology.getNumAtoms()
    print('Number of atoms in cofactor:', ncof)
    print('Call Modeller: include cofactor')
    modeller.add(molcof_ommtopology, molcof_positions)


print('Read ligand 1:')
mollig1 = Molecule.from_file(lig1sdffile, file_format='SDF',
                             allow_undefined_stereo=True)
ligandmolecules.append(mollig1)
lig1_ommtopology = mollig1.to_topology().to_openmm(ensure_unique_atom_names=True)

#assign the residue name to the ligand, assumes one residue
resname_lig1 = "L1"
for residue in lig1_ommtopology.residues():
    residue.name = resname_lig1
pos = mollig1.conformers[0].to('angstrom').magnitude
lig1_positions = [Vec3(pos[i][0], pos[i][1], pos[i][2]) for i in range(pos.shape[0])] * angstrom
nlig1 = lig1_ommtopology.getNumAtoms()
print('Number of atoms in ligand 1:', nlig1)
print('Call Modeller: include ligand 1')
modeller.add(lig1_ommtopology, lig1_positions)

if not rbfe:
    # if ABFE, translate the ligand 1 coordinates into the solvent to calculate the
    # bounding box below
    for i in range(nlig1):
        lig1_positions[i] += displacement
else:
    # RBFE mode:
    # read ligand 2 and place it in the solvent
    print('Read ligand 2:')
    mollig2 = Molecule.from_file(lig2sdffile, file_format='SDF',
                                 allow_undefined_stereo=True)
    ligandmolecules.append(mollig2)
    lig2_ommtopology = mollig2.to_topology().to_openmm(ensure_unique_atom_names=True)
    #assign the residue name, assumes one residue
    resname_lig2 = "L2"
    for residue in lig2_ommtopology.residues():
        residue.name = resname_lig2
    pos = mollig2.conformers[0].to('angstrom').magnitude
    lig2_positions = [Vec3(pos[i][0], pos[i][1], pos[i][2]) for i in range(pos.shape[0])] * angstrom
    nlig2 = lig2_ommtopology.getNumAtoms()
    print('Number of atoms in ligand 2:', nlig2)
    for i in range(nlig2):
        lig2_positions[i] += displacement
    print('Call Modeller: include ligand 2')
    modeller.add(lig2_ommtopology, lig2_positions)
    lig2atom_indexes = [ i for i in range(nrcpt+nlig1,nrcpt+nlig1+nlig2)]
    print("Indexes of ligand 2 (starting from 0):", lig2atom_indexes)


print("Calculating system bounding box:")
if not rbfe:
    bbox = boundingBoxSizes(rcpt_positions + lig1_positions)
else:
    bbox = boundingBoxSizes(rcpt_positions + lig1_positions + lig2_positions)
bboxsizes = [ bbox[i][1]-bbox[i][0] for i in range(3) ]
padding = 5. * 1.0*nanometer
xBoxvec = Vec3((bboxsizes[0]+padding)/nanometer, 0., 0.)*nanometer
yBoxvec = Vec3(0.0, (bboxsizes[1]+padding)/nanometer, 0.)*nanometer
zBoxvec = Vec3(0.0, 0.0, (bboxsizes[2]+padding)/nanometer)*nanometer
print("boxVectors:", (xBoxvec,yBoxvec,zBoxvec ))



#bboxfaces = [ bboxsizes[2]*bboxsizes[1], bboxsizes[2]*bboxsizes[0],  bboxsizes[1]*bboxsizes[0] ]
#print("Areas of faces", bboxfaces)
#smallest_direction = 0
#smallest_area = bboxfaces[0]
#for i in range(3):
#    if bboxfaces[i] < smallest_area:
#        smallest_direction = i
#print("Smallest direction", smallest_direction)


    
############################################
#                                          #
#   SET UP FORCEFIELD FOR LIGANDS          #
#                                          #
############################################

print('\nSet up the combined protein + ligand + water system for simulation')
template_gen = None
if ligandforcefield[0:4] == "gaff":
    from openmmforcefields.generators import GAFFTemplateGenerator
    print('Using GAFFTemplateGenerator function for ligands')
    template_gen = GAFFTemplateGenerator(molecules=ligandmolecules, cache=ffcachefile )
elif ligandforcefield[0:6] == "openff":
    from openmmforcefields.generators import SMIRNOFFTemplateGenerator
    print('Call SMIRNOFFTemplateGenerator function for ligands')
    template_gen = SMIRNOFFTemplateGenerator(molecules=ligandmolecules, forcefield=ligandforcefield, cache=ffcachefile )
elif ligandforcefield[0:8] == "espaloma":
    from openmmforcefields.generators import EspalomaTemplateGenerator
    print('Call EspalomaTemplateGenerator function for ligands')
    template_gen = EspalomaTemplateGenerator(molecules=ligandmolecules, forcefield=ligandforcefield, cache=ffcachefile )
else:
    print('Unknown ligand force field %s' % ligandforcefield)
    sys.exit(1)

# Register the SMIRNOFF template generator
# NOTE: forcefield object was initialized (above)
# for the protein + water (subsystem) force field 
# This step adds support for the ligand force field
forcefield.registerTemplateGenerator(template_gen.generator)

if implsolv is None:
    print("Adding solvent and processing system ...")
    modeller.addSolvent(forcefield, boxVectors = (xBoxvec,yBoxvec,zBoxvec ))
    print("Number of atoms in solvated system:", modeller.topology.getNumAtoms())
    system=forcefield.createSystem(modeller.topology, nonbondedMethod = PME, nonbondedCutoff = 0.9*nanometer,
                                   constraints=HBonds, rigidWater = True, removeCMMotion = False, hydrogenMass = hmass*amu)
else:
    print("Solvent model: %s" % implsolv)
    print("Number of atoms in implicit solvent system:", modeller.topology.getNumAtoms())
    print("Processing system ...")
    system=forcefield.createSystem(modeller.topology, nonbondedMethod = NoCutoff,
                                   constraints=HBonds, rigidWater = True, removeCMMotion = False, hydrogenMass = hmass*amu)

with open(xmloutfile, 'w') as output:
    output.write(XmlSerializer.serialize(system))

if pdboutfile is not None:
    PDBFile.writeFile(modeller.topology, modeller.positions,
                      open(pdboutfile,'w'), keepIds=True)


#find the first atom index of the ligand
lig1_start = modeller.topology.getNumAtoms()
for at in modeller.topology.atoms():
    if at.residue.name == resname_lig1:
        if at.index < lig1_start:
            lig1_start = at.index
if lig1_start == modeller.topology.getNumAtoms():
    raise Exception("Error: could not find ligand %s" % resname_lig1)
    
lig1atom_indexes = [ i for i in range(lig1_start,lig1_start+nlig1)]
print("Indexes of ligand 1 (starting from 0):", lig1atom_indexes)

if rbfe:
    print("Indexes of ligand 2 (starting from 0):", [ i for i in range(lig1_start+nlig1,lig1_start+nlig1+nlig2)])
    
print("Box Vectors:", modeller.topology.getPeriodicBoxVectors())
    
today = datetime.today()
print('\n\nDate and time at end:   ', today)
program_end_timer = time()
print('\nTotal compute time %.3f seconds' %  (program_end_timer-program_start_timer))

