Relative Binding Free Energies of SDC1 Peptide Variants Bound to TIAM1
-----------------------------------------------------------------------

In this tutorial, we will use an automated workflow to calculate the relative binding free energies between single-point mutants of the SDC1 peptide binding to the TIAM1 protein. This system is used as a benchmark for the Alchemical Transfer with Coordinate Swapping (ATS) method in [E. Gallicchio. Relative Binding Free Energy Estimation of Congeneric Ligands and Macromolecular Mutants with the Alchemical Transfer with Coordinate Swapping Method. J. Chem. Inf. Model., 2025.](https://doi.org/10.1021/acs.jcim.5c00207).

The workflow prepares and runs relative binding free energy calculations between pairs of SDC1 mutant variants using a dual-coordinate/dual-topology alchemical algorithm. For general background on the Alchemical Transfer Method (ATM), refer to [Azimi, Khuttan, Wu, Pal, Gallicchio. Relative Binding Free Energy Calculations for Ligands with Diverse Scaffolds with the Alchemical Transfer Method.](https://pubs.acs.org/doi/full/10.1021/acs.jcim.1c01129).

This example uses the `protein_mutation` workflow mode, which is designed for pairwise calculations between macromolecular sequence variants represented as PDB files. For small-molecule ligand RBFE calculations, see the [CDK2 example](../cdk2/README.md) instead.

We assume that the latest release of OpenMM and the latest AToM-OpenMM packages and their dependencies are available in a `conda` environment (see the [README](../../../README.md)), and that the `examples` folder is available under `$HOME/AToM-OpenMM/examples`. Adjust this pathname as needed.

## For the Impatient

Run the workflow:

```
cd $HOME/AToM-OpenMM/examples/RBFE/tiam1
bash ./scripts/setup-atm.sh
cd complexes

#run the jobs on a SLURM cluster
for i in tiam1-* ; do ( cd $i && sbatch ./run.sh ) ; done

#without a SLURM cluster
#for i in tiam1-* ; do ( cd $i && bash ./run.sh > ${i}.log 2>&1 ) ; done
```

Collect the results once the simulations have completed:

```
cd $HOME/AToM-OpenMM/examples/RBFE/tiam1/complexes
for i in tiam1-* ; do ( echo $i && tail -2 ${i}.log ) ; done
```

## Tell me More: I Want to Adapt This to My System

### Input Files

This workflow uses PDB files for both the receptor and the alchemical partner structures. The receptor PDB file is stored in the `receptor/` subdirectory, and the PDB files for each partner (wild-type and mutants) are stored in the `ligands/` subdirectory. All input structures are assumed to be fully prepared, including hydrogen atoms.

The workflow's behavior is primarily controlled by the [`setup-settings.sh`](scripts/setup-settings.sh) file in the `scripts/` folder. In `protein_mutation` mode, this file specifies:

- `workflow_mode=protein_mutation`
- The basename of the receptor PDB file (`tiam1`, corresponding to `receptor/tiam1.pdb`).
- The list of transformations in the form `"MUT1 MUT2 RESID"`, where `MUT1` and `MUT2` are the basenames of the partner PDB files in `ligands/`, and `RESID` is the residue number of the mutation site.

For example, the entry `"sdc1wt sdc1Q3E 3"` in `protein_mutation_transformations` pairs `ligands/sdc1wt.pdb` with `ligands/sdc1Q3E.pdb`, using residue 3 to define the alignment atoms.

### System Preparation

This step sets up the alchemical transfer simulations. Run the automated setup script from the example root directory:

```
cd $HOME/AToM-OpenMM/examples/RBFE/tiam1
bash ./scripts/setup-atm.sh
```

In `protein_mutation` mode, the `setup-atm.sh` script performs the following actions:

- Runs `scripts/find_pair_alignment_atoms.py` on the partner PDB files to generate the pair-specific alignment atom file `ligands/pair_alignments.yaml`.
- Creates simulation directories under `complexes/` for each transformation pair, named following the format `tiam1-sdc1wt-sdc1Q3E`.
- Prepares a `run.sh` launch script in each simulation subdirectory using `scripts/run_template.sh` as a template. `run.sh` works as a SLURM batch script or a regular bash script for interactive use.

### The find_pair_alignment_atoms.py Script

The `scripts/find_pair_alignment_atoms.py` script is the key addition for `protein_mutation` mode. For each transformation pair, it reads the specified residue from both partner PDB files and extracts the 1-based ATOM/HETATM record indices of the backbone atoms in the order CA, N, C. These three atoms define the alignment restraint frame used to keep the orientation of the unbound partner in approximate alignment with the bound partner during the alchemical transfer (see [The Ligand Alignment Restraints](https://www.compmolbiophysbc.org/atom-openmm/atom-system-setup#h.vndnoxipr7qs) for more information).

The results are written to `ligands/pair_alignments.yaml`. Each entry maps the job name (e.g. `tiam1-sdc1wt-sdc1Q3E`) to the alignment atom indices for both partners. This file is generated automatically during setup and does not need to be prepared by hand.

If the specified RESID matches backbone atoms in more than one chain within a partner PDB, the script stops with an error rather than guessing. In that case, use a PDB where the mutation-site residue number is unique across chains.

### Execution and Free Energy Analysis

Submit the alchemical Hamiltonian Replica Exchange simulations to a SLURM cluster, or run them sequentially on a local machine:

```
#run the jobs on a SLURM cluster
for i in tiam1-* ; do ( cd $i && sbatch ./run.sh ) ; done

#without a SLURM cluster
#for i in tiam1-* ; do ( cd $i && bash ./run.sh > ${i}.log 2>&1 ) ; done
```

The `run.sh` launch script calls `scripts/run-atm.py` for each pair. The `run-atm.py` application performs the following tasks:

- It loads default settings from `scripts/defaults.yaml`. This file contains settings common to all pairs, such as the alchemical schedule, the maximum runtime (`WALL_TIME`), and the number of perturbation energy samples to collect (`MAX_SAMPLES`). To keep runtime short for this tutorial, `MAX_SAMPLES` is set to 10 per replica. Set `MAX_SAMPLES` to at least 100 for quantitative work and increase `WALL_TIME` and the SLURM job time limit in `run_template.sh` accordingly. See [The AToM Control File Reference](https://www.compmolbiophysbc.org/atom-openmm/atom-control-file) for a full description of the settings.
- It loads the pair-specific alignment atoms from `ligands/pair_alignments.yaml`.
- It finds an optimal initial position for the unbound partner in the solvent bulk to minimize the simulation box size.
- It builds the simulation system, adding solvent and assigning force field parameters (`amber-14` for the receptor, solvent, and protein partners).
- It runs `rbfe_structprep` to equilibrate and anneal the system.
- It runs `rbfe_production` to collect perturbation energy data using the variable-displacement alchemical Hamiltonian Replica Exchange algorithm.
- It calls `calculate_uwham` to estimate the relative binding free energy using UWHAM thermodynamic reweighting.

After or while production is underway, view the alchemical trajectories with VMD using the provided script. For example:

```
cd $HOME/AToM-OpenMM/examples/RBFE/tiam1/complexes/tiam1-sdc1wt-sdc1Q3E
vmd -f tiam1-sdc1wt-sdc1Q3E_0.pdb `/bin/ls -v r*/*xtc` -e vmd.in
```

### Key Differences from the Ligand RBFE Workflow

Compared to the small-molecule CDK2 example:

- Partner structures are PDB files in `ligands/`, not SDF files.
- Alignment atoms are specific to each transformation pair and are defined by the backbone CA, N, C atoms of the mutation-site residue, rather than being shared across all ligands from a common reference.
- The alchemical variable atoms are the side-chain atoms at the mutation site, not a complete small molecule.
- The `RCPT_CHAIN_NAME` setting in `scripts/defaults.yaml` must match the chain ID of the receptor in your prepared PDB. For this TIAM1 system it is set to `B`. Most standard preparations use chain `A`.

## Credits

Adapted from:  Emilio Gallicchio <emilio.gallicchio@gmail.com>

The algorithms to find the alignment atoms and the optimal displacements are adapted from the [ATM](https://github.com/EricChen521/atm) package by Eric Chen @EricChen521.

The `protein_mutation` workflow mode and the `find_pair_alignment_atoms.py` script were developed to extend the workflow to macromolecular mutant systems.
