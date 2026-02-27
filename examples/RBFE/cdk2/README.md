Relative Binding Free Energies of a set of ligands of the CDK2 Kinase
---------------------------------------------------------------------

In this tutorial, we will use an automated workflow to calculate the relative binding free energies between six congeneric ligands binding to the CDK2 receptor from the paper: [Wang et al. (2013) Modeling Local Structural Rearrangements Using FEP/REST: Application to Relative Binding Affinity Predictions of CDK2 Inhibitors](https://pubs.acs.org/doi/10.1021/ct300911a).

The workflow is set to prepare and run eight relative binding free energy calculations between a selection of pairs of the six CDK2 complexes using a dual-coordinate/dual-topology alchemical algorithm. For more information, refer to [Azimi, Khuttan, Wu, Pal, Gallicchio]. Relative Binding Free Energy Calculations for Ligands with Diverse Scaffolds with the Alchemical Transfer Method.](https://pubs.acs.org/doi/full/10.1021/acs.jcim.1c01129) and [ E. Gallicchio. Relative Binding Free Energy Estimation of Congeneric Ligands and Macromolecular Mutants with the Alchemical Transfer with Coordinate Swapping Method ](https://doi.org/10.1021/acs.jcim.5c00207).

We assume that the latest release of OpenMM and the latest AToM-OpenMM packages and their dependencies are available in a `conda` environment (see the [README](../../../README.md), and that the `examples` folder is available under `$HOME/AToM-OpenMM/examples`. Adjust this pathname as needed.

## For the Impatient

Run the workflow:

```
cd $HOME/AToM-OpenMM/examples/RBFE/cdk2
bash ./scripts/setup-atm.sh
cd complexes

#run the jobs on a SLURM cluster
for i in cdk2-* ; do ( cd $i && sbatch ./run.sh ) ; done

#without a SLURM cluster
#for i in cdk2-* ; do ( cd $i && bash ./run.sh > ${i}.log 2>&1 ) ; done
```

Collect the results once the simulations have completed:

```
cd $HOME/AToM-OpenMM/examples/RBFE/cdk2
for i in cdk2-* ; do ( echo $i && tail -2 ${i}.log ) ; done
```

## Tell me More: I Want to Adapt This to my System

### Input Files

This automated workflow assumes that the `sdf` files of the ligand are stored in the `ligands/` subdirectory and the PDB file of the receptor is stored in the `receptor/` subdirectory. It is assumed that the receptor structure is fully prepared, including hydrogen atoms, etc. The congeneric ligands are assumed to be properly docked to the receptor binding site, and their common rigid substructures are aligned.

The workflow's behavior is primarily controlled by the `setup-settings.sh` file in the `scripts/` folder. This file specifies:

- The basename of the receptor PDB file (`cdk2`, in this case, corresponding to the `cdk2.pdb` file of the receptor in the `receptor/` folder).
- The basename of the reference ligand SDF file (`H1Q`, in this case, corresponding to the `H1Q.sdf` file in the `ligands/` folder).
- The indices (starting from 1) of the three alignment atoms of the reference ligand according to the order in the SDF file above. The alignment atoms are used to keep the orientation of the unbound ligand in approximate alignment with the bound ligand during the alchemial transfer calculation (see [The Ligand Alignment Restraints](https://www.compmolbiophysbc.org/atom-openmm/atom-system-setup#h.vndnoxipr7qs), for more information). The workflow automatically assigns the closest atoms of the other ligands to the alignment atoms of the reference ligand.
- The list of ligand pairs to process.

### System preparation

This step sets up the alchemical transfer simulations. The automated setup script below reads the parameters from the `setup-settings.sh` file.

```
cd $HOME/AToM-OpenMM/examples/RBFE/cdk2
bash ./scripts/setup-atm.sh
```

The `setup-atm.sh` script performs the following actions:

- Runs the `scripts/find_alignment_atoms.py` in the `ligands/` folder that sets the alignment atoms of the ligands.
- Creates simulation directories for the alchemical transfer simulation for each ligand pair under the `complexes/` subdirectory. Simulations directories are named following the format `cdk2-H1Q-H1R`, where `cdk2` is the receptor name and `H1Q` and `H1R` are the ligand names.
- Prepares a `run.sh` launch script in each simulation subdirectory using the `scripts/run_template.sh` as a template. `run.sh` is designed to work as a SLURM batch script or a regular `bash` script for interactive work. Edit `scripts/run_template.sh` to change its behavior.
- The main function of `run.sh` is to execute the `run-atm.py` application (saved under `scripts/') described below.

### Execution and Free Energy Analysis

This step submits the alchemical parallel molecular dynamics simulations to a SLURM cluster. Alternatively, it runs them sequentially on the local machine.  

```
#run the jobs on a SLURM cluster
for i in cdk2-* ; do ( cd $i && sbatch ./run.sh ) ; done

#without a SLURM cluster
#for i in cdk2-* ; do ( cd $i && bash ./run.sh > ${i}.log 2>&1 ) ; done
```

The `run.sh` launch script runs the `scripts/run-atm.py` application, which is the main workhorse of the workflow (run `python scripts/run-atm.py` without arguments to see the arguments it accepts). The `run-atm.py` application performs the following tasks:

- It loads the default options from the `scripts/defaults.yaml` configuration file. This file contains settings (such as the alchemical schedule) common to all ligand pairs. Edit this file to change the default settings. Among other things, `scripts/defaults.yaml` sets the simulation length in terms of the maximum runtime (`WALL_TIME`) and the number of perturbation energy samples to collect (`MAX_SAMPLES`). To minimize runtime, `MAX_SAMPLES` is set to 10 per replica in this tutorial. Set `MAX_SAMPLES` to at least 100 for quantitative work (more is always better, but also more computationally expensive) and increase `WALL_TIME` and the SLURM job max time in `run_template.sh` accordingly. See [The AToM Control File Reference](https://www.compmolbiophysbc.org/atom-openmm/atom-control-file) for a description of the settings in `scripts/defaults.yaml`.
- It loads the ligands' alignment atoms from the `ligands/alignments.yaml` file generated by the `scripts/find_alignment_atoms.py` app above. You can prepare this file in other ways, including by hand, to customize the choice of alignment atoms.
- It finds an optimal initial position for the unbound ligand in the solvent bulk. The algorithm is designed to minimize the simulation box size to optimize performance.
- It builds the simulation system, including adding the solvent and assigning force field parameters (`amber-14` for the receptor and the solvent, and `gaff` for the ligands in this tutorial).
- It sets the simulation control parameters specific to the ligand pair, such as the indices of the ligands and other things. It saves the specific settings together with the default settings in a `.yaml` file in the simulation directory.
- It runs the `rbfe_structprep` application to equilibrate and anneal the system.
- It runs the `rbfe_production` application to run the alchemical Hamiltonian Replica Exchange simulation to collect the perturbation energy data.
- It calls the `calculate_uwham` routine to estimate the relative binding free energy using the UWHAM thermodynamic reweighting method from the collected perturbation energy samples.
