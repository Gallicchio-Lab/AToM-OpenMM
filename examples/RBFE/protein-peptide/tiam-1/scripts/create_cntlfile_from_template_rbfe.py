# Usage: python create_cntlfile_from_template.py <options>
# Emilio Gallicchio, 6/2023

import os, sys, fileinput
import numpy as np
from time import time

# following for argument passing tools
import argparse

# OpenMM components
from openmm import Vec3
from openmm.app import PDBFile, Topology
from openmm.unit import nanometer, angstrom

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
parser.add_argument('--lig1chainname',  required=True,  type=str, default=None,
                    help='Chain id of peptide 1')
parser.add_argument('--lig2chainname',  required=True,  type=str, default=None,
                    help='Chain id of peptide 2')
parser.add_argument('--lig1refatoms',  required=True,  type=str, default=None,
                    help='Reference atoms of ligand 1')
parser.add_argument('--lig2refatoms',  required=True,  type=str, default=None,
                    help='Reference atoms of ligand 2')
parser.add_argument('--vsiteResidues',  required=True,  type=str, default=None,
                    help='Job name')
parser.add_argument('--cntlfileout',  required=True,  type=str, default=None,
                    help='Output control file')
parser.add_argument('--lig1varatoms',  required=True,  type=str, default=None,
                    help='Variable atoms of ligand 1')
parser.add_argument('--lig2varatoms',  required=True,  type=str, default=None,
                    help='Variable atoms of ligand 2')

args = vars(parser.parse_args())

def get_indexes_from_query(topology, query):
    indexes = []
    for atom in topology.atoms():
        if eval(query):
            indexes.append(atom.index)
    indexes.sort()
    return indexes

def cm_from_indexes(topology, positions, indexes):
    cm = Vec3(0,0,0)*nanometer
    n = 0
    for atom in topology.atoms():
        if atom.index in indexes:
            cm += positions[atom.index]
            n += 1
    cm = cm/float(n)
    return cm

def cm_from_query(topology, positions, query):
    indexes = get_indexes_from_query(topology, query)
    return cm_from_indexes(topology, positions, indexes)

def list_to_string(alist):
    s = ""
    if len(alist) < 1:
        return s
    if len(alist) > 1:
        s = str(alist[0])
        for i in alist[1:] :
            s = s + "," + str(i)
        return s
    else:
        return str(alist[0]) + ","

# Pull data from command line into local variables
pdbfile = args['systemPDBFile']
displ = [float(r) for r in args['displacement'].split()]
templatefile = args['templatein']
jobname = args['jobname']
lig1chainname = args['lig1chainname']
lig2chainname = args['lig2chainname']
lig1refatoms = [int(i)-1 for i in args['lig1refatoms'].split(',')]
lig2refatoms = [int(i)-1 for i in args['lig2refatoms'].split(',')]
lig1varatoms = [int(i)-1 for i in args['lig1varatoms'].split(',')]
lig2varatoms = [int(i)-1 for i in args['lig2varatoms'].split(',')]
vsiteresids = [r for r in args['vsiteResidues'].split()]
cntlfile = args['cntlfileout']

#reads the topology
pdb = PDBFile(pdbfile )
positions = pdb.positions
ommtopology = pdb.topology

#list of restrained atoms (Calpha atoms of the receptor)
restr_atoms = []
for atom in ommtopology.atoms():
    if atom.name == "CA" and atom.residue.chain.id == "B":
        restr_atoms.append(atom.index)
print("Positionally restrained atoms:", restr_atoms)
                
#list of ligand atoms
lig1_atoms = []
for at in ommtopology.atoms():
    if at.residue.chain.id == lig1chainname:
        lig1_atoms.append(at.index)
print("Lig1 atoms:", lig1_atoms)

lig2_atoms = []
for at in ommtopology.atoms():
    if at.residue.chain.id == lig2chainname:
        lig2_atoms.append(at.index)
print("Lig2 atoms:", lig2_atoms)

#use the C-alpha atoms of the terminal residues of the peptides as ligand CM
query = 'atom.residue.chain.id == "%s" and atom.residue.id in ["8", "7", "6"] and atom.name == "CA"' % lig1chainname
lig1cmindexes = get_indexes_from_query(ommtopology, query)
lig1cm = cm_from_indexes(ommtopology, positions, lig1cmindexes)
print("Lig1 CM atoms:", lig1cmindexes)
print("Lig1 CM:", lig1cm)
query = 'atom.residue.chain.id == "%s" and atom.residue.id in ["8", "7", "6"] and atom.name == "CA"' % lig2chainname
lig2cmindexes = get_indexes_from_query(ommtopology, query)
lig2cm = cm_from_indexes(ommtopology, positions, lig2cmindexes)
print("Lig2 CM atoms:", lig2cmindexes)
print("Lig2 CM:", lig2cm)
#list of receptor vsite atoms (Calpha atoms of selected residues)
vsite_rcpt_atoms = []
for residue in ommtopology.residues():
    if residue.id in vsiteresids:
        for atom in residue.atoms():
            if atom.name == "CA":
                vsite_rcpt_atoms.append(atom.index)
print("Rcpt Vsite atoms:", vsite_rcpt_atoms)
cmr = cm_from_indexes(ommtopology, positions, vsite_rcpt_atoms)
print("Rcpt CM:",  cmr)

displacement = lig2cm - lig1cm
print("Displacement = ", displacement)
offset = lig1cm - cmr
print("LigOffset = ", offset)

lig1_varatoms = [ lig1_atoms[i] for i in lig1varatoms ]
lig2_varatoms = [ lig2_atoms[i] for i in lig2varatoms ]

lig1_attachatom = lig1refatoms[0] + lig1_atoms[0]
lig2_attachatom = lig2refatoms[0] + lig2_atoms[0]

with open(cntlfile, 'w') as outf:
    for l in fileinput.input(files = templatefile):
        l = l.replace('<JOBNAME>', jobname)
        d = displacement/angstrom
        l = l.replace('<DISPLX>', str(d[0]))
        l = l.replace('<DISPLY>', str(d[1]))
        l = l.replace('<DISPLZ>', str(d[2]))
        d = offset/angstrom
        l = l.replace('<OFFX>', str(d[0]))
        l = l.replace('<OFFY>', str(d[1]))
        l = l.replace('<OFFZ>', str(d[2]))
        l = l.replace('<LIG1ATOMS>', list_to_string(lig1_atoms))
        l = l.replace('<LIG2ATOMS>', list_to_string(lig2_atoms))
        l = l.replace('<LIG1VARATOMS>', list_to_string(lig1_varatoms))
        l = l.replace('<LIG2VARATOMS>', list_to_string(lig2_varatoms))
        l = l.replace('<LIG1ATTACHATOM>', str(lig1_attachatom))
        l = l.replace('<LIG2ATTACHATOM>', str(lig2_attachatom))
        l = l.replace('<LIG1CMATOMS>', list_to_string(lig1cmindexes))
        l = l.replace('<LIG2CMATOMS>', list_to_string(lig2cmindexes))
        l = l.replace('<REFERENCEATOMS1>', list_to_string(lig1refatoms) )
        l = l.replace('<REFERENCEATOMS2>', list_to_string(lig2refatoms) )
        l = l.replace('<VSITERECEPTORATOMS>', list_to_string(vsite_rcpt_atoms))
        l = l.replace('<RESTRAINEDATOMS>', list_to_string(restr_atoms))
        outf.write(l)

with open('vmd.in', 'w') as outf:
    for l in fileinput.input(files = 'vmd_template.in'):
        l = l.replace('<LIG1ATTACHATOM>', str(lig1_attachatom))
        l = l.replace('<LIG2ATTACHATOM>', str(lig2_attachatom))
        outf.write(l)

