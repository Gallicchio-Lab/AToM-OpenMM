#! python

# Usage: python make_atm_rbfe_system_frompdb.py <options>
# Emilio Gallicchio, 5/2023 adapted from code by Bill Swope, 11/2021

############################################
#                                          #
#   IMPORTS                                #
#                                          #
############################################

# following for argument passing tools
import argparse

# following for date and time
from datetime import datetime

# following for timing code components
from time import time

# OpenMM components
from openmm import XmlSerializer
from openmm.app import AmberPrmtopFile, AmberInpcrdFile, PDBFile
from openmm.app import PME, HBonds
from openmm.unit import nanometer, amu


def make_system(prmtopfile, crdfile, xmloutfile, pdboutfile, hmass, nonbondedCutoff, switchDistance, verbose=False):
    #####################################################
    #   Echo out the suppliable input parameters        #
    #####################################################
    print('Generate ATM OpenMM System from Amber files (prmtop and inpcrd)')
    today = datetime.today()
    program_start_timer = time()
    print('\nDate and time at start: ', today.strftime('%c'))

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
    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=nonbondedCutoff*nanometer,
                                                constraints=HBonds, hydrogenMass = float(hmass) * amu, switchDistance=switchDistance*nanometer)

    with open(xmloutfile, 'w') as output:
        output.write(XmlSerializer.serialize(system))

    if pdboutfile is not None:
        PDBFile.writeFile(prmtop.topology, inpcrd.positions, 
                        open(pdboutfile,'w'), keepIds=True )
        
    today = datetime.today()
    print('\n\nDate and time at end:   ', today)
    program_end_timer = time()
    print('\nTotal compute time %.3f seconds' %  (program_end_timer-program_start_timer))
    return system


def main():
    whatItDoes = """
    Produces an .xml file with the OpenMM's system for ATM binding
    free energy calculations from Amber's topology and coordinate files
    produced by tleap.

    Emilio Gallicchio 6/2023
    """
    parser = argparse.ArgumentParser(description=whatItDoes)

    # Required input
    parser.add_argument('--AmberPrmtopinFile', required=True,  type=str, default=None, dest='prmtopfile',
                        help='Amber topology file')
    parser.add_argument('--AmberInpcrdinFile', required=True,  type=str, default=None, dest='crdfile',
                        help='Amber coordinate file')
    parser.add_argument('--systemXMLoutFile',  required=True,  type=str, default=None, dest='xmloutfile',
                        help='Name of the XML file where to save the System')
    parser.add_argument('--systemPDBoutFile', required=True, type=str, default=None, dest='pdboutfile',
                        help='Name of the PDB file where to output to system')

    # Optional input
    parser.add_argument('--hmass', required=False, type=float,
                        default=1.0,
                        help='Hydrogen mass, set it to 1.5 amu to use a 4 fs time-step')
    parser.add_argument('--nonbondedCutoff', required=False, type=float,
                        default=0.9,
                        help='Nonbonded cutoff, default is 0.9 nm')
    parser.add_argument('--switchDistance', required=False, type=float,
                        default=0.0,
                        help='Switch distance, default is 0.0 nm')

    # Arguments that are flags
    parser.add_argument('--verbose', required=False, action='store_true',
                        help='Get more output with this flag')
    args = vars(parser.parse_args())
    make_system(**args)


if __name__ == "__main__":
    main()
