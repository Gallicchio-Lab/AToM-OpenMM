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
                    help='Job name')
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
        
#list of vsite atoms (Calpha atoms of selected residues)
vsite_atoms = []
for residue in ommtopology.residues():
    if residue.id in vsiteresids:
        for atom in residue.atoms():
            if atom.name == "CA":
                vsite_atoms.append(atom.index)

#list of ligand atoms
lig1_atoms = []
for residue in ommtopology.residues():
    if residue.name == lig1resname:
        for atom in residue.atoms():
            lig1_atoms.append(atom.index)

with open(cntlfile, 'w') as outf:
    for l in fileinput.input(files = templatefile):
        l = l.replace('<JOBNAME>', jobname)
        l = l.replace('<DISPLX>', str(displ[0]))
        l = l.replace('<DISPLY>', str(displ[1]))
        l = l.replace('<DISPLZ>', str(displ[2]))
        l = l.replace('<LIGATOMS>', list_to_string(lig1_atoms))
        l = l.replace('<VSITERECEPTORATOMS>', list_to_string(vsite_atoms))
        l = l.replace('<RESTRAINEDATOMS>', list_to_string(restr_atoms))
        outf.write(l)
