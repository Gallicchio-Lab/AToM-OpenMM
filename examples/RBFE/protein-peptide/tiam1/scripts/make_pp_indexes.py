# Thin CLI wrapper around atom_openmm.utils.AtomUtils.make_pp_indexes
# Adapted from Emilio Gallicchio, 7/2025

import argparse

from openmm.app import PDBFile

from atom_openmm.utils.AtomUtils import make_pp_indexes


def list_to_string(alist):
    s = ""
    if len(alist) < 1:
        return s
    if len(alist) > 1:
        s = str(alist[0])
        for i in alist[1:]:
            s = s + "," + str(i)
        return s
    return str(alist[0]) + ","


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--ipdb",
        required=True,
        type=str,
        default=None,
        help="PDB file for topology",
    )
    parser.add_argument(
        "--chainLig1",
        required=True,
        type=str,
        default=None,
        help='Chain id of first ligand, such as "B"',
    )
    parser.add_argument(
        "--residLig1",
        required=True,
        type=str,
        default=None,
        help='Residue id of mutated residue of first ligand, such as "8"',
    )
    parser.add_argument(
        "--chainLig2",
        required=True,
        type=str,
        default=None,
        help='Chain id of second ligand, such as "C"',
    )
    parser.add_argument(
        "--residLig2",
        required=True,
        type=str,
        default=None,
        help='Residue id of mutated residue of second ligand, such as "8"',
    )
    args = vars(parser.parse_args())

    topology = PDBFile(args["ipdb"]).topology
    indexes = make_pp_indexes(
        topology=topology,
        chainLig1=args["chainLig1"],
        chainLig2=args["chainLig2"],
        residLig1=args["residLig1"],
        residLig2=args["residLig2"],
    )

    print("LIGAND1_ATOMS = ", list_to_string(indexes["LIGAND1_ATOMS"]))
    print("LIGAND2_ATOMS = ", list_to_string(indexes["LIGAND2_ATOMS"]))
    print("LIGAND1_VAR_ATOMS = ", list_to_string(indexes["LIGAND1_VAR_ATOMS"]))
    print("LIGAND2_VAR_ATOMS = ", list_to_string(indexes["LIGAND2_VAR_ATOMS"]))
    print("LIGAND1_COMMON_ATOMS = ", list_to_string(indexes["LIGAND1_COMMON_ATOMS"]))
    print("LIGAND2_COMMON_ATOMS = ", list_to_string(indexes["LIGAND2_COMMON_ATOMS"]))
    print("LIGAND1_ATTACH_ATOM = ", indexes["LIGAND1_ATTACH_ATOM"])
    print("LIGAND2_ATTACH_ATOM = ", indexes["LIGAND2_ATTACH_ATOM"])
    print(
        "ALIGN_LIGAND1_REF_ATOMS = ",
        list_to_string(indexes["ALIGN_LIGAND1_REF_ATOMS"]),
    )
    print(
        "ALIGN_LIGAND2_REF_ATOMS = ",
        list_to_string(indexes["ALIGN_LIGAND2_REF_ATOMS"]),
    )
    print("LIGAND1_CM_ATOMS =", list_to_string(indexes["LIGAND1_CM_ATOMS"]))
    print("LIGAND2_CM_ATOMS =", list_to_string(indexes["LIGAND2_CM_ATOMS"]))


if __name__ == "__main__":
    main()
