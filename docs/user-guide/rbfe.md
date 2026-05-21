# Relative Binding Free Energy

Relative binding free energy (RBFE) workflows estimate free energy differences between two related bound states. In AToM-OpenMM, the newer RBFE examples use a YAML-based setup that prepares paired alchemical systems, runs alchemical Hamiltonian replica exchange, and analyzes the two alchemical legs with UWHAM.

The most useful local examples to study are:

- `examples/RBFE/cdk2`: small-molecule ligand pairs with generated alignment atoms.
- `examples/RBFE/protein-peptide/tiam1`: protein-peptide variants where a mutation residue defines the common and variable atoms.

The CDK2 and protein-peptide workflows are the best templates for the current `setup-settings.sh` plus `defaults.yaml` pattern.

## When to Use RBFE

Use RBFE when you care about a difference between two states: ligand A versus ligand B, peptide variant A versus peptide variant B, or one molecular variant versus another. RBFE is usually the natural choice for congeneric ligand series, local chemical changes, or mutation studies where a perturbation network can be planned.

Use ABFE instead when each ligand should be treated independently or when there is no meaningful relative perturbation path.

## Small-Molecule Ligand Pairs

The CDK2 workflow is the current small-molecule RBFE template. It assumes:

- A prepared receptor PDB in `receptor/`.
- Ligand SDF files in `ligands/`.
- Docked or modeled ligand poses in the binding site.
- A reference ligand and three reference alignment atoms.
- A list of ligand pairs to process.

The key setup file is `scripts/setup-settings.sh`. For CDK2, it defines the receptor basename, the reference ligand, the reference ligand alignment atoms, and the ligand pairs:

```text
receptor=cdk2
ref_ligand=H1Q
ref_ligand_alignment_atoms="14,21,18"
ligands=( "H1Q H1R" "H1Q H1S" ... )
```

When `scripts/setup-atm.sh` runs, it first calls `scripts/find_alignment_atoms.py` to generate `ligands/alignments.yaml`. The generated alignment file lets `run-atm.py` assign `ALIGN_LIGAND1_REF_ATOMS` and `ALIGN_LIGAND2_REF_ATOMS` for each ligand pair. These atoms are used by the orientation and roll components of the ligand alignment restraint.

The CDK2 defaults use variable displacement: `ALIGN_KF_SEP` is zero, while `ALIGN_K_THETA` and `ALIGN_K_PSI` restrain relative orientation. During production, the workflow swaps the positions of the two ligands based on their current anchor-atom separation, confines the bound ligand to the binding-site region, and uses a short-range exclusion potential to keep the displaced ligand from interacting with the receptor.

## Protein-Peptide and Mutation Workflows

For peptide variants or mutation-like transformations, study `examples/RBFE/protein-peptide/tiam1`. This workflow uses the same overall launch structure as CDK2, but each entry in `scripts/setup-settings.sh` includes a mutation residue:

```text
ligands=( "sdc1E4Q sdc1wt 4"
          "sdc1wt  sdc1A8F 8" )
```

The two peptide or protein-partner PDB files are passed to `run-atm.py` as `LIG1` and `LIG2`. The `--mutationResid` value is written as `MUTATION_RESID` and used by the shared `make_pp_indexes()` helper to identify:

- The full atom lists for the two partners.
- The common and variable atoms around the mutation site.
- The attachment atoms.
- The alignment reference atoms derived from the mutation-residue backbone.

By default, the helper treats `N`, `CA`, `C`, `O`, and `H` as backbone atoms. Set `PP_BACKBONE` in `scripts/defaults.yaml` if the partner naming or chemistry requires a different backbone definition.

Use this protein-peptide pattern when the perturbation is better described by a residue or local molecular variant than by small-molecule SDF alignment. The important requirement is that the receptor and partner structures are prepared consistently and that the mutation residue can be identified in both partners.

## Running an RBFE Workflow

For the CDK2 small-molecule workflow:

```bash
cd $HOME/AToM-OpenMM/examples/RBFE/cdk2
bash ./scripts/setup-atm.sh
cd complexes

# On a SLURM cluster
for i in cdk2-*; do ( cd "$i" && sbatch ./run.sh ); done

# Or without SLURM
# for i in cdk2-*; do ( cd "$i" && bash ./run.sh > "${i}.log" 2>&1 ); done
```

For the TIAM1 protein-peptide workflow:

```bash
cd $HOME/AToM-OpenMM/examples/RBFE/protein-peptide/tiam1
bash ./scripts/setup-atm.sh
cd complexes

# On a SLURM cluster
for i in tiam1-*; do ( cd "$i" && sbatch ./run.sh ); done

# Or without SLURM
# for i in tiam1-*; do ( cd "$i" && bash ./run.sh > "${i}.log" 2>&1 ); done
```

Each generated `run.sh` calls `scripts/run-atm.py`, which loads `scripts/defaults.yaml`, adds pair-specific atom selections and displacement information, writes `<jobname>.yaml`, runs `rbfe_structprep`, runs `rbfe_production`, and performs UWHAM analysis when production data are available.

## Outputs to Check

Inside each `complexes/<jobname>/` directory, the main outputs are:

| Output | Use |
| --- | --- |
| `<jobname>.yaml` | Final merged options for this pair. |
| `<jobname>.pdb`, `<jobname>_sys.xml` | Prepared dual-topology or paired system. |
| `<jobname>_0.xml`, `<jobname>_0.pdb` | Prepared state used to start production. |
| `r*/<jobname>.out` | Perturbation-energy samples from each replica. |
| `<jobname>.log` | Runtime log and final relative free-energy estimate. |
| `<jobname>.png` | UWHAM quality-control plot when requested. |
| `vmd.in` | VMD helper file generated from the template when available. |

The final log lines report `DG` in kcal/mol, its estimated uncertainty, the two leg free energies, and the number of samples used after discarding initial samples.

## Planning and Adapting

For small-molecule RBFE, plan the perturbation network before setup. Choose ligand pairs that are scientifically meaningful and that have reliable bound poses. The CDK2 workflow can generate alignment atoms from a reference ligand, but you should still inspect the chosen atoms and the resulting structures. Poor alignment atoms can lead to unstable restraints or unhelpful ligand orientations.

For mutation or protein-peptide RBFE, make sure the partner structures use consistent residue numbering and atom naming. The mutation residue should identify the local change cleanly in both partners. If the default backbone definition does not match the input structures, set `PP_BACKBONE` explicitly.

For either RBFE type, tune these files first:

| File | What to adjust |
| --- | --- |
| `scripts/setup-settings.sh` | Receptor name, ligand pairs, peptide pairs, mutation residues, and reference alignment settings. |
| `scripts/defaults.yaml` | Alchemical schedule, force field, receptor chains, restraints, displacement behavior, and production length. |
| `scripts/run_template.sh` | Conda environment activation, SLURM resources, time limit, and command-line options passed to `run-atm.py`. |

Tutorial defaults are intentionally short. For production work, increase `MAX_SAMPLES`, increase `WALL_TIME`, keep scheduler limits consistent, and compare quality-control plots and log summaries across the perturbation network before interpreting final rankings.
