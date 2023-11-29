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
from openmm import XmlSerializer
from openmm.app import AmberPrmtopFile, AmberInpcrdFile, PDBFile
from openmm.app import PME, HBonds

from openmm.unit import Quantity
from openmm.unit import angstrom, nanometer, nanometers, picoseconds, amu
from openmm import unit

from sys import stdout

# System calls - to invoke sdfTagTool from the python code
import os, sys

print('Generate ATM OpenMM System from Amber files (prmtop and inpcrd)')
today = datetime.today()
print('\nDate and time at start: ', today.strftime('%c'))


whatItDoes = """
Produces an .xml file with the OpenMM's system for ATM binding
free energy calculations from Amber's topology and coordinate files
produced by tleap.

Emilio Gallicchio 6/2023
"""

#####################################################
#                                                   #
#   PASS COMMAND LINE ARGUMENTS TO LOCAL VARIABLES  #
#                                                   #
#####################################################


program_start_timer = time()
parser = argparse.ArgumentParser(description=whatItDoes)

# Required input
parser.add_argument('--AmberPrmtopinFile', required=True,  type=str, default=None,
                    help='Amber topology file')
parser.add_argument('--AmberInpcrdinFile', required=True,  type=str, default=None,
                    help='Amber coordinate file')
parser.add_argument('--systemXMLoutFile',  required=True,  type=str, default=None,
                    help='Name of the XML file where to save the System')
parser.add_argument('--systemPDBoutFile', required=True, type=str, default=None,
                    help='Name of the PDB file where to output to system')

# Optional input
parser.add_argument('--hmass', required=False, type=float,
                    default=1.0,
                    help='Hydrogen mass, set it to 1.5 amu to use a 4 fs time-step')

# Arguments that are flags
parser.add_argument('--verbose', required=False, action='store_true',
                    help='Get more output with this flag')

args = vars(parser.parse_args())

# Pull data from command line into local variables
prmtopfile = args['AmberPrmtopinFile']
crdfile = args['AmberInpcrdinFile']
xmloutfile = args['systemXMLoutFile']
pdboutfile = args['systemPDBoutFile']
hmass = float(args['hmass'])

# flag for printing (verbose) 
flagverbose = args['verbose']

#####################################################
#   Echo out the suppliable input parameters        #
#####################################################


print('\nUser-supplied input parameters')
print('Amber prmtop file:               ', prmtopfile)
print('Amber coordinate file:           ', crdfile)
print('Topology PDB output file:        ', pdboutfile)
print('System XML output file:          ', xmloutfile)
print('Hydrogen mass:                   ', hmass)  


############################################
#                                          #
#   READ THE AMBER SYSTEM                  #
#                                          #
############################################

prmtop = AmberPrmtopFile(prmtopfile)
inpcrd = AmberInpcrdFile(crdfile)
system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.9*nanometer,
                                               constraints=HBonds, hydrogenMass = hmass * amu)

with open(xmloutfile, 'w') as output:
    output.write(XmlSerializer.serialize(system))

if pdboutfile is not None:
    PDBFile.writeFile(prmtop.topology, inpcrd.positions, 
                      open(pdboutfile,'w'))
    
today = datetime.today()
print('\n\nDate and time at end:   ', today)
program_end_timer = time()
print('\nTotal compute time %.3f seconds' %  (program_end_timer-program_start_timer))

