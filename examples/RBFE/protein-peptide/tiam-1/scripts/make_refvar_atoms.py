# Usage: python create_refatoms.py <options>
# Emilio Gallicchio, 6/2023

import os, sys, fileinput, re
import numpy as np
from time import time

# following for argument passing tools
import argparse

# OpenMM components
from openmm.app import PDBFile, Topology

def list_to_string(alist):
    s = str(alist[0])
    for i in alist[1:] :
        s = s + "," + str(i)
    return s

def refatoms(top, mut):
    pattern = re.compile('[A-Z]([1-9][0-9]?)[A-Z]p?')
    matches = pattern.search(mut)
    resnum = int(matches[1])
    ref = [-1,-1,-1]
    for atom in top.atoms():
        rid = -1
        if int(atom.residue.id) == resnum:
            if atom.name == "CA":
                rid = 0
            elif atom.name == "N":
                rid = 1
            elif atom.name == "C":
                rid = 2
        if rid >= 0:
            ref[rid] = atom.index+1
    return ref

def varatoms(top, mut):
    pattern = re.compile('[A-Z]([1-9][0-9]?)[A-Z]')
    matches = pattern.search(mut)
    resnum = int(matches[1])
    var = []
    backbone = """["N","CA","C","O", "OXT", "H","H2", "H3","HA", "DN", "DCA", "DC", "DO", "LPOA", "LPOB",  "DOT1", "DOT2", "LPT1", "LPT2", "LPT3", "LPT4"]"""
    for atom in top.atoms():
        if int(atom.residue.id) == resnum:
            if not atom.name in backbone:
                var.append(atom.index+1)
    return var

parser = argparse.ArgumentParser()

# Required input
parser.add_argument('--mutations', required=True, type=str, default=None,
                    help='List of mutations "wt mut1 E55A, wt mut2 V51A, ... " ')
args = vars(parser.parse_args())

# Pull data from command line into local variables
mutinput = args['mutations'].split(',')
lig1 = []
lig2 = []
mutations = []
for m in mutinput:
    (l1, l2, mut) = m.split()
    lig1.append(l1)
    lig2.append(l2)
    mutations.append(mut)

print("ligands=(")
for i in range(len(mutations)):
    l1 = lig1[i]
    l2 = lig2[i]
    mut = mutations[i]
    print('\"' + l1 + ' ' + l2 + '\"')
print(")")
    
refwt = []
refmut = []
varwt = []
varmut = []
for i in range(len(mutations)):
    l1 = lig1[i]
    l2 = lig2[i]
    mut = mutations[i]
    pdbfilewt = l1 + ".pdb"
    topwt = PDBFile( pdbfilewt ).topology
    pdbfilemut = l2 + ".pdb"
    topmut =  PDBFile( pdbfilemut ).topology
    refwt.append(refatoms(topwt, mut).copy())
    refmut.append(refatoms(topmut, mut).copy())
    varwt.append(varatoms(topwt, mut).copy())
    varmut.append(varatoms(topmut, mut).copy())

print("ref_atoms=(")
for i in range(len(mutations)):
    print('\"' + list_to_string(refwt[i]) + ' ' + list_to_string(refmut[i]) + '\"')
print(")")

print("declare -A variable_region=(")
for i in range(len(mutations)):
    l1 = lig1[i]
    l2 = lig2[i]
    print('[\"' + l1 + ' ' + l2 + '\"]=' + '\"' + list_to_string(varwt[i]) + '    ' + list_to_string(varmut[i]) + '\"')
print(")")
