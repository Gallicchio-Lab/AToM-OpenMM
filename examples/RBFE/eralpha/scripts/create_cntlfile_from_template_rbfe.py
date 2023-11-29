# Usage: python create_cntlfile_from_template.py <options>
# Emilio Gallicchio, 6/2023

import os, sys, fileinput
import numpy as np
from time import time

# following for argument passing tools
import argparse

# OpenMM components
from openmm.app import PDBFile, Topology

parser = argparse.ArgumentParser()

# Required input
parser.add_argument('--systemPDBFile', required=True, type=str, default=None,
                    help='Name of the PDB file with the system')
parser.add_argument('--displacement',  required=True,  type=str, default=None,
                    help='string with displacement vector in angstroms like "22. 0.0 0.0" ')
parser.add_argument('--templatein',  required=True,  type=str, default=None,
                    help='Template cntl file')
parser.add_argument('--jobname',  required=True,  type=str, default=None,
                    help='Job name')
parser.add_argument('--lig1resname',  required=True,  type=str, default=None,
                    help='Residue name of ligand 1')
parser.add_argument('--lig2resname',  required=True,  type=str, default=None,
                    help='Residue name of ligand 1')
parser.add_argument('--lig1refatoms',  required=True,  type=str, default=None,
                    help='Reference atoms of ligand 1')
parser.add_argument('--lig2refatoms',  required=True,  type=str, default=None,
                    help='Reference atoms of ligand 2')
parser.add_argument('--vsiteResidues',  required=True,  type=str, default=None,
                    help='Job name')
parser.add_argument('--cntlfileout',  required=True,  type=str, default=None,
                    help='Output control file')

args = vars(parser.parse_args())

def list_to_string(alist):
    s = str(alist[0])
    for i in alist[1:] :
        s = s + ", " + str(i)
    return s

# Pull data from command line into local variables
pdbfile = args['systemPDBFile']
displ = [float(r) for r in args['displacement'].split()]
templatefile = args['templatein']
jobname = args['jobname']
lig1resname = args['lig1resname']
lig2resname = args['lig2resname']
lig1refatoms = [int(i)-1 for i in args['lig1refatoms'].split(',')]
lig2refatoms = [int(i)-1 for i in args['lig2refatoms'].split(',')]
vsiteresids = [r for r in args['vsiteResidues'].split()]
cntlfile = args['cntlfileout']

#reads the topology
pdb = PDBFile(pdbfile )
positions = pdb.positions
ommtopology = pdb.topology

#list of restrained atoms (Calpha atoms)
restr_atoms = []
for atom in ommtopology.atoms():
    if atom.name == "CA":
        restr_atoms.append(atom.index)
print("Positionally restrained atoms:", restr_atoms)
        
#list of vsite atoms (Calpha atoms of selected residues)
vsite_atoms = []
for residue in ommtopology.residues():
    if residue.id in vsiteresids:
        for atom in residue.atoms():
            if atom.name == "CA":
                vsite_atoms.append(atom.index)
print("VSite atoms:", vsite_atoms)
                
#list of ligand atoms
lig1_atoms = []
for residue in ommtopology.residues():
    if residue.name == lig1resname:
        for atom in residue.atoms():
            lig1_atoms.append(atom.index)
        break #reads the first match only
print("Lig1 atoms:", lig1_atoms)

lig2_atoms = []
for residue in ommtopology.residues():
    if residue.name == lig2resname:
        for atom in residue.atoms():
            #sometimes lig1 and lig2 have the same residue name
            if not atom.index in lig1_atoms: 
                lig2_atoms.append(atom.index)
print("Lig2 atoms:", lig2_atoms)

#use the first reference atom as the CM
lig1cm = lig1_atoms[lig1refatoms[0]]
lig2cm = lig2_atoms[lig2refatoms[0]]

with open(cntlfile, 'w') as outf:
    for l in fileinput.input(files = templatefile):
        l = l.replace('<JOBNAME>', jobname)
        l = l.replace('<DISPLX>', str(displ[0]))
        l = l.replace('<DISPLY>', str(displ[1]))
        l = l.replace('<DISPLZ>', str(displ[2]))
        l = l.replace('<LIG1ATOMS>', list_to_string(lig1_atoms))
        l = l.replace('<LIG2ATOMS>', list_to_string(lig2_atoms))
        l = l.replace('<LIG1CMATOMS>', str(lig1cm) + ",")
        l = l.replace('<LIG2CMATOMS>', str(lig2cm) + ",")
        l = l.replace('<REFERENCEATOMS1>', list_to_string(lig1refatoms) )
        l = l.replace('<REFERENCEATOMS2>', list_to_string(lig2refatoms) )
        l = l.replace('<VSITERECEPTORATOMS>', list_to_string(vsite_atoms))
        l = l.replace('<RESTRAINEDATOMS>', list_to_string(restr_atoms))
        outf.write(l)
