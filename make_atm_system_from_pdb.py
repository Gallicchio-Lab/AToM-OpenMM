# Usage: python make_atm_rbfe_system_frompdb.py <options>
# Emilio Gallicchio, 5/2023 adapted from code by Bill Swope, 11/2021

############################################
#                                          #
#   IMPORTS                                #
#                                          #
############################################

# following for argument passing tools
import argparse

# following used to generate a subdirectory named after the ligand
# for the output files
import subprocess

# following for date and time
from datetime import datetime

# following for timing code components
from time import time

# OpenMM components
# from simtk.openmm.app import Modeller, ForceField, Simulation
# from simtk.openmm.app import PDBReporter, StateDataReporter, PDBFile
# from simtk.openmm.app import PME, HBonds
from openmm import XmlSerializer
from openmm.app import Modeller, ForceField, Simulation
from openmm.app import PDBReporter, StateDataReporter, PDBFile
from openmm.app import PME, HBonds

# from simtk.openmm import Platform, MonteCarloBarostat, LangevinMiddleIntegrator
# from simtk.openmm import Vec3
from openmm import Platform, MonteCarloBarostat, LangevinMiddleIntegrator
from openmm import Vec3

# from simtk.unit import Quantity, bar, kelvin
# from simtk.unit import angstrom, nanometer, nanometers, picoseconds
# from simtk import unit
from openmm.unit import Quantity
from openmm.unit import angstrom, nanometer, nanometers, picoseconds, amu
from openmm import unit

from sys import stdout

# System calls - to invoke sdfTagTool from the python code
import os, sys

# OpenFF and OpenMM components for ligand force field parameters
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit.topology import Molecule

# OpenFF components for SystemGenerator (for input entire ssytem in pdb file)
from openmm import app
from openmmforcefields.generators import SystemGenerator


############################################
#                                          #
#   ACTUAL CODE                            #
#                                          #
############################################

print('Generate ATM RBFE OpenMM System')
today = datetime.today()
print('\nDate and time at start: ', today.strftime('%c'))


whatItDoes = """
Produces an .xml file with the OpenMM's system for ATM relative binding
free energy calcuations.

The user supplies a pdb file that contains protein, cofactors, ligands, 
water, ions, probably prepared by a different software package 
(Maestro, OpenEye).  Cofactors and ligands have to be described in an
sdf file.  Then force fields are assigned to all components

Emilio Gallicchio 5/2023
adapted from simProteinLigandWater.py script by Bill Swope 11/2021
"""

#####################################################
#                                                   #
#   PASS COMMAND LINE ARGUMENTS TO LOCAL VARIABLES  #
#                                                   #
#####################################################


program_start_timer = time()
parser = argparse.ArgumentParser(description=whatItDoes)

# Required input
parser.add_argument('--systemPDBinFile', required=True,  type=str, default=None,
                    help='Prepared system in PDB format')
parser.add_argument('--ligandsSDFFile',  required=True,  type=str, default=None,
                    help='SDF file with all ligands and cofactors')
parser.add_argument('--LIG1resid',  required=True,  type=int, default=None,
                    help='Residue id of the first ligand')
parser.add_argument('--systemXMLoutFile',  required=True,  type=str, default=None,
                    help='Name of the XML file where to save the System')
parser.add_argument('--systemPDBoutFile', required=True, type=str, default=None,
                    help='Name of the PDB file where to output to system')

# Optional input
parser.add_argument('--LIG1refatoms',  required=False,  type=str, default=None,
                    help='string of reference atoms of the first ligand "1 2 3" ')
parser.add_argument('--LIG2resid',  required=False,  type=int, default=None,
                    help='Residue id of the second ligand')
parser.add_argument('--LIG2refatoms',  required=False,  type=str, default=None,
                    help='string of reference atoms of the first ligand "1 2 3"')
parser.add_argument('--proteinForceField', required=False, type=str,
                    default='amber14-all.xml',
                    help='Force field for protein` amber14-all.xml ')
parser.add_argument('--solventForceField', required=False, type=str,
                    default='amber14/tip3p.xml',
                    help='Force field for solvent/ions` amber14/tip3p.xml ')
parser.add_argument('--ligandForceField', required=False, type=str,
                    default='openff-2.0.0',
                    help='Force field for ligand:  openff-2.0.0')
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

#catch abfe or rbfe
rbfe = False
if args['LIG2resid']  is not None:
    if ( (args['LIG1refatoms'] is not None)   and
         (args['LIG2resid'] is not None)      and
         (args['LIG2refatoms'] is not None) ):
        rbfe = True
    else:
        print("Error: missing information for RBFE setup")
        print(parser.print_help())
        sys.exit(1)
    
# Pull data from command line into local variables
systempdbfile = args['systemPDBinFile']
ligandsdffile = args['ligandsSDFFile']
xmloutfile = args['systemXMLoutFile']
pdboutfile = args['systemPDBoutFile']
lig1resid = args['LIG1resid']
if rbfe:
    lig2resid = args['LIG2resid']
    lig1refatoms = [int(i) for i in args['LIG1refatoms'].split()]
    lig2refatoms = [int(i) for i in args['LIG2refatoms'].split()]

# force fields
proteinforcefield = args['proteinForceField']
solventforcefield = args['solventForceField']
ligandforcefield = args['ligandForceField']
ffcachefile = args['forcefieldJSONCachefile']

hmass = float(args['hmass'])

# flag for printing (verbose) 
flagverbose = args['verbose']

#####################################################
#   Echo out the suppliable input parameters        #
#####################################################


print('\nUser-supplied input parameters')
print('System PDB file name:               ', systempdbfile)
print('Protein force field:                ', proteinforcefield)
print('Solvent/ion force field             ', solventforcefield )
print('Ligands file name:                  ', ligandsdffile)
print('Ligand force field:                 ', ligandforcefield)
if rbfe:
    print('Residue id of first ligand:         ', lig1resid)
    print('Reference atoms of first ligand:    ', lig1refatoms)
    print('Residue id of second ligand:        ', lig2resid)
    print('Reference atoms of second ligand    ', lig2refatoms)
else:
    print('Residue id of ligand:               ', lig1resid)
print('Topology PDB output file:           ', pdboutfile)
print('System XML output file:             ', xmloutfile)
print('Force field cache file:             ', ffcachefile)


print('\nAvailable small molecule OpenFF force fields for ligand:')
print(SMIRNOFFTemplateGenerator.INSTALLED_FORCEFIELDS)



############################################
#                                          #
#   READ AND CHARACTERIZE SYSTEM           #
#                                          #
############################################

print('\nSystem characteristics')
print('File name:      ', systempdbfile)
systempdb = PDBFile(systempdbfile)
print('Number of atoms in pdb file:', len(systempdb.positions))
print('Box vectors:', systempdb.topology.getPeriodicBoxVectors())

############################################
#                                          #
#   READ AND CHARACTERIZE LIGANDS           #
#                                          #
############################################

lig1atomindexes = []
for atom in systempdb.topology.atoms():
    if int(atom.residue.id) == lig1resid:
        lig1atomindexes.append(atom.index)

if rbfe:
    reflig1atomindexes = [lig1atomindexes[0]+i-1 for i in lig1refatoms]
    lig2atomindexes = []
    for atom in systempdb.topology.atoms():
        if int(atom.residue.id) == lig2resid:
            lig2atomindexes.append(atom.index)
    reflig2atomindexes = [lig2atomindexes[0]+i-1 for i in lig2refatoms]
        
    print("Atom indexes of first ligand (starting from 0):")
    print(lig1atomindexes)
    print("Atom indexes of the reference atoms of first ligand:")
    print(reflig1atomindexes)

    print("Atom indexes of second ligand (starting from 0):")
    print(lig2atomindexes)
    print("Atom indexes of the reference atoms of second ligand:")
    print(reflig2atomindexes)

    #obtain the displacement from the distance between the two first reference atoms
    displacement = systempdb.positions[reflig2atomindexes[0]] - systempdb.positions[reflig1atomindexes[0]]
    print("ATM Displacement vector:", displacement)
else:
    print("Atom indexes of ligand (starting from 0):")
    print(lig1atomindexes)

print('\nLigands characteristics')
print('Ligand file name:           ', ligandsdffile)
ligmolecules = Molecule.from_file(ligandsdffile, file_format='SDF',
                                  allow_undefined_stereo=True)
print(ligmolecules)

periodic_forcefield_kwargs = {'nonbondedMethod': PME, 'nonbondedCutoff': 0.9*nanometer}
forcefield_kwargs={'constraints' : HBonds, 'rigidWater' : True,
                   'removeCMMotion' : False, 'hydrogenMass' : hmass*amu }

system_generator = SystemGenerator(forcefields=[ proteinforcefield, solventforcefield],
                                   small_molecule_forcefield=ligandforcefield,
                                   forcefield_kwargs=forcefield_kwargs,
                                   periodic_forcefield_kwargs=periodic_forcefield_kwargs,
                                   cache=ffcachefile)
print(system_generator.forcefield.getUnmatchedResidues(systempdb.topology))
print(system_generator.forcefield.generateTemplatesForUnmatchedResidues(systempdb.topology))
system=system_generator.create_system(systempdb.topology,
                                      molecules=ligmolecules)

with open(xmloutfile, 'w') as output:
    output.write(XmlSerializer.serialize(system))

if pdboutfile is not None:
    PDBFile.writeFile(systempdb.topology, systempdb.positions, 
                      open(pdboutfile,'w'))

today = datetime.today()
print('\n\nDate and time at end:   ', today)
program_end_timer = time()
print('\nTotal compute time %.3f seconds' %  (program_end_timer-program_start_timer))

