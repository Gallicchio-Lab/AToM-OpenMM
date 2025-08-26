# AToM-OpenMM Protein-Peptide RBFE Workflow Tutorial 

This workflow tutorial sets up and runs the relative binding free energies between the Syndecan-1 derived peptide mutants and the TIAM-1 PDZ domain described in the paper: [Relative Binding Free Energy Estimation of Congeneric Ligands and Macromolecular Mutants with the Alchemical Transfer Method with Coordinate Swapping](https://pubs.acs.org/doi/10.1021/acs.jcim.5c00207).

The setup used here is sligtly different than in the paper. Here the binding site region is defined in terms of the distance between the centroid of the C-alpha atoms of the first three residues of the peptides (0, 1, and 2 following the numbering in the paper) and the centroid of the C-alpha atoms of the binding site residues of the TIAM-1 PDZ receptor (see below). 


## OpenMM/AToM-OpenMM Requirements

This tutorial requires OpenMM 8.4 or newer. More specifically a version of OpenMM built after the commit 3645e1e ([Variable distance-based displacements for ATMForce](https://github.com/openmm/openmm/commit/3645e1eeb7cea37bc7af1291e0fd672a85d49db9)). It also requires the latest version of the AToM-OpenMM package.

Follow the installation instructions in the main [README](https://github.com/Gallicchio-Lab/AToM-OpenMM/blob/master/README.md). Here we assume that all required packages are installed under a miniforge3 environment activated using: 
```
. ${HOME}/miniforge3/bin/activate
```
Modify the command to suit your specific environment.

## Setup

### Short story

Assuming the AToM-OpenMM repo is available under your home directory, do:
```
cd ${HOME}/AToM-OpenMM/examples/RBFE/protein-peptide/tiam-1
. ${HOME}/miniforge3/bin/activate
bash scripts/setup-atm.sh
```
This creates a `complexes` folder with subdirectories each corresponding to the RBFE simulation for a pair of mutants. See below for execution and analysis.

### Longer story

The PDB of the TIAM-1 PDZ receptor is stored under the `receptor/` directory and PDBs of the peptide mutants are stored under the `ligands/` directory. We created these files using Maestro starting with the structure of the wild-type complex (4GVD) and applying the mutations.

The `scripts` folder hold the necessary scripts and template input files. In particular, The `setup-atm.sh` script above is a bash script that orchestrates the setup of each RBFE simulation with the help of a couple of python scripts.

The first, called `make_refvar_atoms.py`, analyzes the peptides to find the atom indexes of the sidechain being mutated. For example:
```
python ${scripts_dir}/make_refvar_atoms.py --mutations "sdc1wt sdc1A0F A8F, sdc1A0F sdc1wt F8A"
```
generates settings for the mutation of alanine 8 (0 in the paper) to phenylalanine and of the reverse mutations, where `sdc1wt.pdb` and `sdc1A0F.pdb` are the corresponding PDB files of the variants with alanine and phenylalanine, respectively.

The script calls AToM-OpenMM's `make_atm_system_from_rcpt_lig` utility to create the simulation system for each mutation by displacing the second peptide by the prescribed displacement vector (the `displacement` setting in the `setup-atm.sh` script) and solvating with water and ions.

Next, The `scripts/create_cntlfile_from_template_rbfe.py` script generates AToM-OpenMM's control file (ending in `.cntl`) for each simulation. It adds restrains for the C-alpha atoms of the receptor with a tolerance of 5.0 Angstroms and it determines the atom indexes of the set of atoms that define the centroids for the binding site definition among other things.

## Execution

The setup generates `slurm` runfiles for each prepared simulations under `complexes/`. Do:
```
cd complexes
for i in $(/bin/ls -d */) ; do ( dir=${i%%/} && cd $dir && sbatch run.sh ) ; done
```
to launch the simulations. Replace `run.sh` with `run4.sh` to use 4 GPUs for each job.

Each job is set to run for 16.6 hours (1,000 minutes). Use `squeue` to monitor the jobs. When the jobs have completed, the command above can be used to relaunch them to obtain more data.

### Get the Relative Free Energies

After the jobs collected a reasonable amount of data, run the analysis script in the complexes directory to obtain the relative free energies:
```
cd complexes
bash ./free_energies.sh
```








