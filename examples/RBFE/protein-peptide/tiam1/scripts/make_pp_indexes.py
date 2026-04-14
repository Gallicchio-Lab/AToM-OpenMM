# Usage:
#python ~/Dropbox/src/workflows-scripts/atm-protein-protein-rbfe-from-CHARMM/scripts/make_pp_indexes.py --ipdb tiam1-sdc1wt-sdc1A0V.pdb --chainLig1 B --chainLig2 C --residLig1 8 --residLig2 8 --rcptCMquery 'atom.residue.chain.id == "A" and atom.name == "CA" and atom.residue.id in ["850", "856", "857", "858", "859", "860", "861", "862", "912", "915", "916", "920"]' --PosRestrquery 'atom.name == "CA" and not (atom.residue.chain.id == "C") and  not ( (atom.residue.chain.id == "B") and (atom.residue.id == "7" or atom.residue.id == "8"))'
# Emilio Gallicchio, 7/2025

import os, sys, fileinput, re
import numpy as np
from time import time

# following for argument passing tools
import argparse

# OpenMM components
from openmm.app import PDBFile, Topology

def get_indexes_from_top(topology, query):
    indexes = []
    for atom in topology.atoms():
        if eval(query):
            indexes.append(atom.index)
    indexes.sort()
    return indexes

def get_indexes_from_residue(residue, query):
    indexes = []
    for atom in residue.atoms():
        if eval(query):
            indexes.append(atom.index)
    indexes.sort()
    return indexes

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

def commonatoms(topology,  lig1_atoms, varatoms1, lig2_atoms, varatoms2):
    ca1 = sorted(list(set(lig1_atoms)-set(varatoms1)))
    ca2 = sorted(list(set(lig2_atoms)-set(varatoms2)))
    commonatoms1 = []
    commonatoms2 = []
    atoms = list(topology.atoms())

    for i in ca1:
        #print(f"ca1: {i}, name: {atoms[i].name} residue: {atoms[i].residue.id} chain {atoms[i].residue.chain.id}")
        #search corresponding atom in second ligand
        found = False
        for j in ca2:
            #print(f"ca2: {j}, name: {atoms[j].name} residue: {atoms[j].residue.id} chain {atoms[j].residue.chain.id}")
            if atoms[i].residue.id == atoms[j].residue.id and atoms[i].name == atoms[j].name:
                commonatoms1.append(i)
                commonatoms2.append(j)
                found = True
                break
        if not found:
            print(f"commonatoms(): error: did not find corresponding atom for atom {i}, name: {atoms[i].name} residue: {atoms[i].residue.id} chain {atoms[i].residue.chain.id}")
            return None, None
    if len(ca1) != len(commonatoms1) or len(ca2) != len(commonatoms2):
        print("commonatoms(): error common atoms do not match")
        return None, None
    return commonatoms1, commonatoms2

def make_pp_indexes(topology,
                    chainLig1,
                    chainLig2,
                    residLig1,
                    residLig2,
                    rcpt_cm_atoms = None,
                    pos_restrained_atoms = None):

    backbone = """["N","NT", "CA","CAY", "CAT", "C", "CY", "O", "OY", "OXT", "H", "H1", "H2", "H3", "HT1", "HT2", "HT3", "HNT", "HY1", "HY2", "HY3", "HA", "DN", "DCA", "DC", "DO", "LPOA", "LPOB",  "DOT1", "DOT2", "LPT1", "LPT2", "LPT3", "LPT4"]"""

    #LIGAND1_ATOMS
    chainid = chainLig1
    lig1chain = None
    for chain in topology.chains():
        if chain.id == chainid:
            lig1chain = chain
    ligand1_atoms = sorted([ atom.index for atom in lig1chain.atoms() ])

    #LIGAND2_ATOMS
    chainid = chainLig2
    lig2chain = None
    for chain in topology.chains():
        if chain.id == chainid:
            lig2chain = chain
    ligand2_atoms = sorted([ atom.index for atom in lig2chain.atoms() ])

    #LIGAND1_VAR_ATOMS
    resid = residLig1
    res1 = None
    for residue in lig1chain.residues():
        if residue.id == resid:
            res1 = residue
    ligand1_var_atoms = get_indexes_from_residue(res1, f"not atom.name in {backbone}")

    #LIGAND2_VAR_ATOMS
    resid = residLig2
    res2 = None
    for residue in lig2chain.residues():
        if residue.id == resid:
            res2 = residue
    ligand2_var_atoms = get_indexes_from_residue(res2, f"not atom.name in {backbone}")

    #common atoms
    ligand1_common_atoms, ligand2_common_atoms = commonatoms(topology, ligand1_atoms, ligand1_var_atoms, ligand2_atoms, ligand2_var_atoms)
    if not ligand1_common_atoms:
        raise ValueError("Error in commonatoms()")

    #reference atoms of ligand1
    chainid = chainLig1
    resid = residLig1
    aname = 'CA'
    query = f""" atom.residue.chain.id == '{chainid}' and atom.residue.id == '{resid}' and atom.name == '{aname}' """
    i = get_indexes_from_top(topology, query)[0]
    ligand1_attachment_atom = i
    aname = 'N'
    query = f""" atom.residue.chain.id == '{chainid}' and atom.residue.id == '{resid}' and atom.name == '{aname}' """
    j = get_indexes_from_top(topology, query)[0]
    aname = 'C'
    query = f""" atom.residue.chain.id == '{chainid}' and atom.residue.id == '{resid}' and atom.name == '{aname}' """
    k = get_indexes_from_top(topology, query)[0]
    ligand1_ref_atoms = [i-ligand1_atoms[0], j-ligand1_atoms[0], k-ligand1_atoms[0]]

    #reference atoms of ligand2
    chainid = chainLig2
    resid = residLig2
    aname = 'CA'
    query = f""" atom.residue.chain.id == '{chainid}' and atom.residue.id == '{resid}' and atom.name == '{aname}' """
    i = get_indexes_from_top(topology, query)[0]
    ligand2_attachment_atom = i
    aname = 'N'
    query = f""" atom.residue.chain.id == '{chainid}' and atom.residue.id == '{resid}' and atom.name == '{aname}' """
    j = get_indexes_from_top(topology, query)[0]
    aname = 'C'
    query = f""" atom.residue.chain.id == '{chainid}' and atom.residue.id == '{resid}' and atom.name == '{aname}' """
    k = get_indexes_from_top(topology, query)[0]
    ligand2_ref_atoms = [i-ligand2_atoms[0], j-ligand2_atoms[0], k-ligand2_atoms[0]]

    indexes = {}
    indexes["LIGAND1_ATOMS"] = ligand1_atoms
    indexes["LIGAND2_ATOMS"] = ligand2_atoms
    indexes["LIGAND1_VAR_ATOMS"] = ligand1_var_atoms
    indexes["LIGAND2_VAR_ATOMS"] = ligand2_var_atoms
    indexes["LIGAND1_COMMON_ATOMS"] = ligand1_common_atoms
    indexes["LIGAND2_COMMON_ATOMS"] = ligand2_common_atoms
    indexes["LIGAND1_ATTACH_ATOM"] = ligand1_attachment_atom
    indexes["LIGAND2_ATTACH_ATOM"] = ligand2_attachment_atom
    indexes["ALIGN_LIGAND1_REF_ATOMS"] = ligand1_ref_atoms
    indexes["ALIGN_LIGAND2_REF_ATOMS"] = ligand2_ref_atoms
    indexes["LIGAND1_CM_ATOMS"] = [ ligand1_attachment_atom ]
    indexes["LIGAND2_CM_ATOMS"] = [ ligand2_attachment_atom ]
    if rcpt_cm_atoms is not None:
        indexes["RCPT_CM_ATOMS"] = rcpt_cm_atoms
    if pos_restrained_atoms is not None:
        indexes["POS_RESTRAINED_ATOMS"] = pos_restrained_atoms

    return indexes

def main():
    parser = argparse.ArgumentParser()

    # Required input
    parser.add_argument('--ipdb', required=True, type=str, default=None,
                        help='PDB file for topology')
    parser.add_argument('--chainLig1', required=True, type=str, default=None,
                        help='Chain id of first ligand, such as "B"')
    parser.add_argument('--residLig1', required=True, type=str, default=None,
                        help='Residue id of mutated residue of first ligand, such as "8"')
    parser.add_argument('--chainLig2', required=True, type=str, default=None,
                        help='Chain id of second ligand, such as "C"')
    parser.add_argument('--residLig2', required=True, type=str, default=None,
                        help='Residue id of mutated residue of second ligand, such as "8"')
    parser.add_argument('--PosRestrquery', required=False, type=str, default='atom.name == "CA"',
                        help='Query to locate Calpha atoms to restrain, such as: atom.name == "CA"')
    parser.add_argument('--rcptCMquery', required=False, type=str, default=None,
                        help='Query to locate Calpha atoms to restrain, such as: atom.name == "CA"')
    args = vars(parser.parse_args())

    topology = PDBFile(args['ipdb']).topology

    rcpt_cm_atoms = None
    if args['rcptCMquery'] is not None:
        rcpt_cm_atoms = get_indexes_from_top(topology, args['rcptCMquery'])

    pos_restrained_atoms = None
    if args['PosRestrquery'] is not None:
        pos_restrained_atoms = get_indexes_from_top(topology, args['PosRestrquery'])

    indexes = make_pp_indexes(topology = topology,
                              chainLig1 = args['chainLig1'],
                              chainLig2 = args['chainLig2'],
                              residLig1 = args['residLig1'],
                              residLig2 = args['residLig2'],
                              rcpt_cm_atoms = rcpt_cm_atoms,
                              pos_restrained_atoms = pos_restrained_atoms)

    print("LIGAND1_ATOMS = ", list_to_string(indexes["LIGAND1_ATOMS"]))
    print("LIGAND2_ATOMS = ", list_to_string(indexes["LIGAND2_ATOMS"]))
    print("LIGAND1_VAR_ATOMS = ", list_to_string(indexes["LIGAND1_VAR_ATOMS"]))
    print("LIGAND2_VAR_ATOMS = ", list_to_string(indexes["LIGAND2_VAR_ATOMS"]))
    print("LIGAND1_COMMON_ATOMS = ", list_to_string(indexes["LIGAND1_COMMON_ATOMS"]))
    print("LIGAND2_COMMON_ATOMS = ", list_to_string(indexes["LIGAND2_COMMON_ATOMS"]))
    print("LIGAND1_ATTACH_ATOM = ", indexes["LIGAND1_ATTACH_ATOM"])
    print("LIGAND2_ATTACH_ATOM = ", indexes["LIGAND2_ATTACH_ATOM"])
    print("ALIGN_LIGAND1_REF_ATOMS = ", list_to_string(indexes["ALIGN_LIGAND1_REF_ATOMS"]))
    print("ALIGN_LIGAND2_REF_ATOMS = ", list_to_string(indexes["ALIGN_LIGAND2_REF_ATOMS"]))
    print("LIGAND1_CM_ATOMS =", list_to_string(indexes["LIGAND1_CM_ATOMS"]))
    print("LIGAND2_CM_ATOMS =", list_to_string(indexes["LIGAND2_CM_ATOMS"]))
    if 'RCPT_CM_ATOMS' in indexes:
        print("RCPT_CM_ATOMS =", list_to_string(indexes["RCPT_CM_ATOMS"]))
    if 'POS_RESTRAINED_ATOMS' in indexes:
        print("POS_RESTRAINED_ATOMS =", list_to_string(indexes["POS_RESTRAINED_ATOMS"]))

if __name__ == "__main__":
    main()
