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

from openmm import XmlSerializer
from openmm.app import PDBFile
from openmm.app import PME, HBonds
from openmm.unit import nanometer, amu

# System calls - to invoke sdfTagTool from the python code
import sys

# OpenFF and OpenMM components for ligand force field parameters
from openmmforcefields.generators import SMIRNOFFTemplateGenerator
from openff.toolkit.topology import Molecule

# OpenFF components for SystemGenerator (for input entire ssytem in pdb file)
from openmmforcefields.generators import SystemGenerator



def make_system(
        systempdbfile,
        ligandsdffile,
        lig1resid,
        xmloutfile,
        pdboutfile,
        lig1refatoms=None,
        lig2resid=None,
        lig2refatoms=None,
        proteinforcefield='amber14-all.xml',
        solventforcefield='amber14/tip3p.xml',
        ligandforcefield='openff-2.0.0',
        ffcachefile=None,
        hmass=1.0,
        flagverbose=False
    ):
    print('Generate ATM RBFE OpenMM System')
    today = datetime.today()
    print('\nDate and time at start: ', today.strftime('%c'))
    program_start_timer = time()

    #catch abfe or rbfe
    rbfe = False
    if lig2resid  is not None:
        if ( (lig1refatoms is not None)   and
            (lig2resid is not None)      and
            (lig2refatoms is not None) ):
            rbfe = True
        else:
            print("Error: missing information for RBFE setup. Check the help for the required arguments.")
            sys.exit(1)
        
    if rbfe:
        lig2resid = lig2resid
        if isinstance(lig1refatoms, str):
            lig1refatoms = [int(i) for i in lig1refatoms.split()]
        if isinstance(lig2refatoms, str):
            lig2refatoms = [int(i) for i in lig2refatoms.split()]

    # force fields
    hmass = float(hmass)

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
                        open(pdboutfile,'w'), keepIds=True)

    today = datetime.today()
    print('\n\nDate and time at end:   ', today)
    program_end_timer = time()
    print('\nTotal compute time %.3f seconds' %  (program_end_timer-program_start_timer))


def main():
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
    parser = argparse.ArgumentParser(description=whatItDoes)

    # Required input
    parser.add_argument('--systemPDBinFile', required=True,  type=str, default=None, dest='systempdbfile',
                        help='Prepared system in PDB format')
    parser.add_argument('--ligandsSDFFile',  required=True,  type=str, default=None, dest='ligandsdffile',
                        help='SDF file with all ligands and cofactors')
    parser.add_argument('--LIG1resid',  required=True,  type=int, default=None, dest='lig1resid',
                        help='Residue id of the first ligand')
    parser.add_argument('--systemXMLoutFile',  required=True,  type=str, default=None, dest='xmloutfile',
                        help='Name of the XML file where to save the System')
    parser.add_argument('--systemPDBoutFile', required=True, type=str, default=None, dest='pdboutfile',
                        help='Name of the PDB file where to output to system')

    # Optional input
    parser.add_argument('--LIG1refatoms',  required=False,  type=str, default=None, dest='lig1refatoms',
                        help='string of reference atoms of the first ligand "1 2 3" ')
    parser.add_argument('--LIG2resid',  required=False,  type=int, default=None, dest='lig2resid',
                        help='Residue id of the second ligand')
    parser.add_argument('--LIG2refatoms',  required=False,  type=str, default=None, dest='lig2refatoms',
                        help='string of reference atoms of the first ligand "1 2 3"')
    parser.add_argument('--proteinForceField', required=False, type=str, dest='proteinforcefield',
                        default='amber14-all.xml',
                        help='Force field for protein` amber14-all.xml ')
    parser.add_argument('--solventForceField', required=False, type=str, dest='solventforcefield',
                        default='amber14/tip3p.xml',
                        help='Force field for solvent/ions` amber14/tip3p.xml ')
    parser.add_argument('--ligandForceField', required=False, type=str, dest='ligandforcefield',
                        default='openff-2.0.0',
                        help='Force field for ligand:  openff-2.0.0')
    parser.add_argument('--forcefieldJSONCachefile', required=False, type=str, dest='ffcachefile',
                        default=None,
                        help='Force field ligand cache database')
    parser.add_argument('--hmass', required=False, type=float,
                        default=1.0,
                        help='Hydrogen mass, set it to 1.5 amu to use a 4 fs time-step')

    # Arguments that are flags
    parser.add_argument('--verbose', required=False, action='store_true', dest='flagverbose',
                        help='Get more output with this flag')
    args = vars(parser.parse_args())
    make_system(**args)


if __name__ == "__main__":
    main()
