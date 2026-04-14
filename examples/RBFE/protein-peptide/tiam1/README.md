Relative Binding Free Energies of a set of Syndecan-1 Derived Peptide Variants Binding to the TIAM1 PDZ Domain
---------------------------------------------------------------------------------------------------------------

In this tutorial, we will use an automated workflow to calculate the relative binding free energies between Syndecan-1-derived peptide variants binding to the TIAM1 PDZ receptor from the paper: [Relative Binding Free Energy Estimation of Congeneric Ligands and Macromolecular Mutants with the Alchemical Transfer with Coordinate Swapping Method](https://doi.org/10.1021/acs.jcim.5c00207).

The workflow is set up to prepare and run two relative binding free energy calculations between pairs of Sdc1 peptide variants using a dual-coordinate/dual-topology alchemical algorithm. For more information, refer to [Azimi, Khuttan, Wu, Pal, Gallicchio. Relative Binding Free Energy Calculations for Ligands with Diverse Scaffolds with the Alchemical Transfer Method.](https://pubs.acs.org/doi/full/10.1021/acs.jcim.1c01129) and [ E. Gallicchio. Relative Binding Free Energy Estimation of Congeneric Ligands and Macromolecular Mutants with the Alchemical Transfer with Coordinate Swapping Method ](https://doi.org/10.1021/acs.jcim.5c00207).

We assume that the latest release of OpenMM and the latest AToM-OpenMM packages and their dependencies are available in a `conda` environment (see the [README](../../../README.md), and that the `examples` folder is available under `$HOME/AToM-OpenMM/examples`. Adjust this pathname as needed.

## For the Impatient

Run the workflow:

```
cd $HOME/AToM-OpenMM/examples/RBFE/protein-peptide/tiam1-protein-peptide
bash ./scripts/setup-atm.sh
cd complexes

#run the jobs on a SLURM cluster
for i in tiam1-* ; do ( cd $i && sbatch ./run.sh ) ; done

#without a SLURM cluster
#for i in tiam1-* ; do ( cd $i && bash ./run.sh > ${i}.log 2>&1 ) ; done
```

Collect the results once the simulations have completed:

```
cd $HOME/AToM-OpenMM/examples/RBFE/protein-peptide/tiam1-protein-peptide/complexes
for i in tiam1-* ; do ( echo $i && tail -2 $i/${i}.log ) ; done
```

## Tell me More: I Want to Adapt This to my System

### Input Files

This automated workflow assumes that the PDB files of the peptide variants are stored in the `ligands/` subdirectory and the PDB file of the receptor is stored in the `receptor/` subdirectory. It is assumed that the receptor structure and the peptide structures are fully prepared, including hydrogen atoms, and that the peptide variants share the same residue numbering outside the mutation site.

The workflow's behavior is primarily controlled by the [`setup-settings.sh`](scripts/setup-settings.sh) file in the `scripts/` folder. This file specifies:

- The basename of the receptor PDB file (`tiam1`, in this case, corresponding to the `tiam1.pdb` file of the receptor in the `receptor/` folder).
- The list of peptide pairs to process.
- The residue id of the shared mutation site for each peptide pair. Each entry in `ligands=( ... )` has the form `"PEPTIDE1 PEPTIDE2 MUTATION_RESID"`.

### System preparation

This step sets up the alchemical transfer simulations. The automated setup script below reads the parameters from the `setup-settings.sh` file.

```
cd $HOME/AToM-OpenMM/examples/RBFE/protein-peptide/tiam1-protein-peptide
bash ./scripts/setup-atm.sh
```

The `setup-atm.sh` script performs the following actions:

- Creates simulation directories for the alchemical transfer simulation for each peptide pair under the `complexes/` subdirectory. Simulation directories are named following the format `tiam1-sdc1wt-sdc1E4Q`, where `tiam1` is the receptor name and `sdc1wt` and `sdc1E4Q` are the peptide names.
- Prepares a `run.sh` launch script in each simulation subdirectory using the `scripts/run_template.sh` as a template. `run.sh` is designed to work as a SLURM batch script or a regular `bash` script for interactive work. Edit `scripts/run_template.sh` to change its behavior.
- Copies the `scripts/vmd_template.in` file to each simulation directory so a pair-specific VMD input file can be prepared during execution.
- The main function of `run.sh` is to execute the `run-atm.py` application (saved under `scripts/`) described below.

### Execution and Free Energy Analysis

This step submits the alchemical parallel molecular dynamics simulations to a SLURM cluster. Alternatively, it runs them sequentially on the local machine.  

```
#run the jobs on a SLURM cluster
for i in tiam1-* ; do ( cd $i && sbatch ./run.sh ) ; done

#without a SLURM cluster
#for i in tiam1-* ; do ( cd $i && bash ./run.sh > ${i}.log 2>&1 ) ; done
```

The `run.sh` launch script is launched for each peptide pair. `run.sh` runs the `scripts/run-atm.py` application, which is the main workhorse of the workflow (run `python scripts/run-atm.py` without arguments to see the arguments it accepts). The `run-atm.py` application performs the following tasks:

- It loads the default options from the `scripts/defaults.yaml` configuration file. This file contains settings (such as the alchemical schedule) common to all peptide pairs. Edit this file to change the default settings. Among other things, `scripts/defaults.yaml` sets the simulation length in terms of the maximum runtime (`WALL_TIME`) and the number of perturbation energy samples to collect (`MAX_SAMPLES`). To minimize runtime, `MAX_SAMPLES` is set to 10 per replica in this tutorial. Set `MAX_SAMPLES` to at least 100 for quantitative work (more is always better, but also more computationally expensive) and increase `WALL_TIME` and the SLURM job max time in `run_template.sh` accordingly. See [The AToM Control File Reference](https://www.compmolbiophysbc.org/atom-openmm/atom-control-file) for a description of the settings in `scripts/defaults.yaml`.
- It finds an optimal initial position for the second peptide in the solvent bulk. The algorithm is designed to minimize the simulation box size to optimize performance.
- It builds the simulation system, including adding the solvent and assigning force field parameters.
- It constructs the receptor internal frame from the generated system PDB and uses that as the source of truth for receptor center-of-mass atoms and positional restraint atoms.
- It calls `make_pp_indexes.py` on the generated system topology to determine the peptide-specific ATM keys, including the peptide atom lists, the common-variable atom split, the attachment atoms, and `ALIGN_LIGAND1_REF_ATOMS` / `ALIGN_LIGAND2_REF_ATOMS` from the mutation-residue `CA`, `N`, and `C` atoms.
- It sets the simulation control parameters specific to the peptide pair and saves the specific settings, along with the default settings, in a `.yaml` file in the simulation directory.
- It runs the `rbfe_structprep` application to equilibrate and anneal the system.
- It runs the `rbfe_production` application to run the alchemical Hamiltonian Replica Exchange simulation to collect the perturbation energy data. The calculations in this tutorial use the variable-displacement algorithm, which (alchemically) swaps the positions of the two peptides based on the current distance between the anchor atoms of the peptides. Additionally, the bound peptide is confined to a spherical binding-site region defined with respect to an internal coordinate frame of the receptor. This allows the complex and the unbound peptide to move freely and independently in the simulation box. A custom short-ranged repulsion potential prevents the unbound peptide from interacting with the receptor.
- It calls the `calculate_uwham` routine to estimate the relative binding free energy using the UWHAM thermodynamic reweighting method from the collected perturbation energy samples.


After or while the production is underway, optionally view the alchemical trajectories with VMD using the provided script. For example:
```
cd $HOME/AToM-OpenMM/examples/RBFE/protein-peptide/tiam1-protein-peptide/complexes/tiam1-sdc1wt-sdc1E4Q
vmd -f tiam1-sdc1wt-sdc1E4Q_0.pdb `/bin/ls -v r*/*xtc` -e vmd.in
```


## Credits

Adapted from Emilio Gallicchio <emilio.gallicchio@gmail.com>

The algorithms to find the optimal displacements are adapted from the [ATM](https://github.com/EricChen521/atm) package by Eric Chen @EricChen521
