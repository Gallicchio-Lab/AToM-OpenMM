# User Guide

The user guide explains how the AToM-OpenMM example workflows are organized and how to adapt them to new systems. The tutorials show complete worked examples; this section focuses on the reusable pieces behind those examples.

The current local workflows use a YAML-based pattern. A small setup file lists the systems to prepare, a shared `defaults.yaml` file holds common simulation settings, and a generated `run.sh` script launches one calculation at a time.

## Workflow Map

Most current examples follow the same sequence:

1. Prepare input structures for the receptor and ligands, peptides, or molecular variants.
2. Edit `scripts/setup-settings.sh` to define the receptor name and the ligands, pairs, or mutation sites to process.
3. Edit `scripts/defaults.yaml` for the shared alchemical schedule, force-field choices, restraint settings, and runtime limits.
4. Run `scripts/setup-atm.sh` to create one calculation directory under `complexes/` for each system or perturbation.
5. Launch each generated `run.sh` script with `sbatch` on a SLURM cluster or with `bash` for local execution.
6. Inspect the generated YAML file, log file, UWHAM estimate, and quality-control plot in each calculation directory.

During execution, `run.sh` calls `scripts/run-atm.py`. That script loads `defaults.yaml`, adds system-specific information such as atom indexes and displacement vectors, writes the merged per-system YAML file, runs structure preparation, runs production sampling, and performs UWHAM analysis when enough samples are available.

## Which Workflow to Use

Use an ABFE workflow when you want the binding free energy of one ligand or molecular partner against a receptor. The newer ABFE example treats the physical ligand as `L1` and a one-atom noninteracting ghost particle as `L2`, so it can reuse the same RBFE-style preparation and production machinery.

Use an RBFE workflow when you want a free energy difference between two related states. The current examples cover small-molecule ligand pairs and protein-peptide variants. In both cases, the workflow builds a paired alchemical system, determines the atom selections needed by ATM, and estimates the relative binding free energy from the two alchemical legs.

## Where to Make Changes

| To change... | Edit... |
| --- | --- |
| Which receptor, ligands, peptide pairs, or mutations are prepared | `scripts/setup-settings.sh` |
| Shared simulation settings and alchemical parameters | `scripts/defaults.yaml` |
| Scheduler options, conda environment activation, or launch command layout | `scripts/run_template.sh` before running setup |
| One generated calculation launch script | `complexes/<system>/run.sh` |
| The exact options used by one prepared calculation | `complexes/<system>/<system>.yaml` |

For routine use, prefer changing `setup-settings.sh`, `defaults.yaml`, and `run_template.sh`, then regenerate calculation directories. The per-system YAML files are most useful as records of what was actually prepared or as debugging targets for one calculation.

## Guide Pages

- [Configuration files](configuration.md): YAML defaults, generated per-system files, and the main runtime keywords.
- [Absolute binding free energy](abfe.md): ABFE-specific setup notes and output interpretation.
- [Relative binding free energy](rbfe.md): RBFE-specific planning, alignment, and perturbation-network notes.

Tutorial examples are configured for short demonstration runs. For quantitative calculations, increase sampling settings such as `MAX_SAMPLES`, increase `WALL_TIME`, and keep the scheduler time limit in `run_template.sh` consistent with those choices.
