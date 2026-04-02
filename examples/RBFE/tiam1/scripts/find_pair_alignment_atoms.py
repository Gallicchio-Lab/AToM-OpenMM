"""
find_pair_alignment_atoms.py -- Extract backbone alignment atom indices for protein mutation pairs.

For each transformation pair specified as "MUT1 MUT2 RESID", this script reads the
corresponding PDB files from the ligands/ directory and extracts the 1-based ATOM/HETATM
record indices of the backbone atoms CA, N, and C at the given mutation-site residue.
These three atoms define the alignment restraint frame used to keep the orientation of
the unbound partner in approximate alignment with the bound partner during the alchemical
transfer calculation.

The results are written to a YAML file (ligands/pair_alignments.yaml) that maps each
job name to the pair-specific alignment atom indices for both partners. This file is
consumed automatically by run-atm.py and does not need to be prepared by hand.

If the specified RESID matches backbone atoms in more than one chain within a partner
PDB, the script raises an error rather than guessing. In that case, use a PDB where the
mutation-site residue number is unique across chains.
"""
import argparse
from pathlib import Path

import yaml


ATOM_ORDER = ("CA", "N", "C")


def _maybe_replace_altloc(existing_altloc, new_altloc):
    preferred = {"": 0, "A": 1}
    return preferred.get(new_altloc, 2) < preferred.get(existing_altloc, 2)


def _collect_backbone_atom_ids(pdb_path, residue_id):
    atoms_by_chain = {}
    atom_index = 0

    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue

            atom_index += 1
            atom_name = line[12:16].strip()
            altloc = line[16].strip()
            chain = line[21].strip()
            resid = line[22:26].strip()

            if resid != residue_id:
                continue

            if atom_name not in ATOM_ORDER:
                continue

            atoms = atoms_by_chain.setdefault(chain, {})
            existing = atoms.get(atom_name)
            if existing is None or _maybe_replace_altloc(existing["altloc"], altloc):
                atoms[atom_name] = {"index": atom_index, "altloc": altloc}

    matching_chains = []
    missing_by_chain = {}
    for chain, atoms in atoms_by_chain.items():
        missing = [atom_name for atom_name in ATOM_ORDER if atom_name not in atoms]
        if missing:
            missing_by_chain[chain] = missing
            continue
        matching_chains.append(chain)

    if not matching_chains:
        missing_summary = ", ".join(
            f"{chain or '[blank]'} missing {missing}"
            for chain, missing in sorted(missing_by_chain.items())
        )
        raise ValueError(
            f"{pdb_path}: residue {residue_id} is missing backbone atoms {ATOM_ORDER}. "
            f"Observed matches: {missing_summary or 'none'}"
        )

    if len(matching_chains) > 1:
        raise ValueError(
            f"{pdb_path}: residue {residue_id} is ambiguous across chains {matching_chains}. "
            "Use a PDB where the mutant partner residue numbering is unique."
        )

    chain = matching_chains[0]
    atoms = atoms_by_chain[chain]
    return [atoms[atom_name]["index"] for atom_name in ATOM_ORDER]


def _parse_pair(pair_spec):
    parts = pair_spec.split()
    if len(parts) != 3:
        raise ValueError(
            f'Invalid protein_mutation pair "{pair_spec}". Expected "MUT1 MUT2 RESID".'
        )

    lig1_name, lig2_name, residue_spec = parts
    residue_id = residue_spec.strip()
    if not residue_id:
        raise ValueError(
            f'Invalid residue spec "{residue_spec}" in pair "{pair_spec}". Expected RESID.'
        )

    return lig1_name, lig2_name, residue_id, residue_spec


def build_pair_alignments(job_prefix, structures_dir, pair_specs, file_extension):
    structures_dir = Path(structures_dir)
    results = {}

    for pair_spec in pair_specs:
        lig1_name, lig2_name, residue_id, residue_spec = _parse_pair(pair_spec)
        jobname = f"{job_prefix}-{lig1_name}-{lig2_name}"
        lig1_path = structures_dir / f"{lig1_name}.{file_extension}"
        lig2_path = structures_dir / f"{lig2_name}.{file_extension}"

        if not lig1_path.exists():
            raise FileNotFoundError(f"Missing first mutant structure: {lig1_path}")
        if not lig2_path.exists():
            raise FileNotFoundError(f"Missing second mutant structure: {lig2_path}")

        lig1_align_atom_ids = _collect_backbone_atom_ids(lig1_path, residue_id)
        lig2_align_atom_ids = _collect_backbone_atom_ids(lig2_path, residue_id)

        results[jobname] = {
            "lig1_name": lig1_name,
            "lig2_name": lig2_name,
            "residue": residue_spec,
            "atom_order": list(ATOM_ORDER),
            "lig1_align_atom_ids": lig1_align_atom_ids,
            "lig2_align_atom_ids": lig2_align_atom_ids,
        }

    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--jobPrefix", required=True, help="Job prefix, usually the receptor basename.")
    parser.add_argument(
        "--structuresDir",
        required=True,
        help="Directory containing the mutant PDB structures.",
    )
    parser.add_argument(
        "--pairAlignmentsYAMLoutFile",
        required=True,
        help="YAML file where pair-specific alignment atoms are written.",
    )
    parser.add_argument(
        "--pair",
        action="append",
        dest="pairs",
        default=[],
        help='Protein mutation pair spec in the form "MUT1 MUT2 RESID".',
    )
    parser.add_argument(
        "--fileExtension",
        default="pdb",
        help="Structure file extension for mutant inputs, default: pdb",
    )
    args = parser.parse_args()

    if not args.pairs:
        parser.error("At least one --pair argument is required.")

    alignments = build_pair_alignments(
        job_prefix=args.jobPrefix,
        structures_dir=args.structuresDir,
        pair_specs=args.pairs,
        file_extension=args.fileExtension.lstrip("."),
    )

    with open(args.pairAlignmentsYAMLoutFile, "w", encoding="utf-8") as handle:
        yaml.dump(alignments, handle, default_flow_style=None, width=1000000, sort_keys=False)
