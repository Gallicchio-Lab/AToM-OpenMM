Relative Binding Free Energies of a set of Antagonists of the Estrogen Receptor α Nuclear Receptor
--------------------------------------------------------------------------------------------------

In this tutorial we will calculate the relative binding free energies between four protein-ligand complexes of the ERα nuclear receptor with four antagonists from the paper: [Azimi, Khuttan, Wu, Pal, Gallicchio. Relative Binding Free Energy Calculations for Ligands with Diverse Scaffolds with the Alchemical Transfer Method.](https://pubs.acs.org/doi/10.1021/acs.jcim.1c01129) This tutorial uses an automated workflow to prepare and run the six relative binding free energy calculations between all pairs of complexes.

It is highly recommended to go through less complex [tutorials](https://github.com/Gallicchio-Lab/AToM-OpenMM/tree/master/examples)before attempting this one. The following assumes familiarity with the terms and procedures introduced in the [TEMOA G1 ABFE tutorial](https://github.com/Gallicchio-Lab/AToM-OpenMM/tree/master/examples/ABFE/temoa-g1) and others.

We assume that the OpenMM and AToM-OpenMM packages and their dependencies are available in a `conda` environment (see the [README](../../../README.md), and that the `examples` folder is available under `$HOME/AToM-OpenMM/examples`. Adjust this pathname as needed.

### System preparation

This tutorial assumes that the [Ambertools](http://ambermd.org) executables are in the search path. They can be conveniently installed under the same conda environment:

Setup the simulation input files. The automated setup script below reads the parameters from the `setup-settings.sh` file.
```
cd $HOME/AToM-OpenMM/examples/RBFE/eralpha
bash ./scripts/setup-atm.sh
```
The automated script runs `antechamber` and `tleap` from the AmberTools package to prepare the ligands and the receptor and to solvated them in a solution box. The `make_atm_system_from_amber` command is used to convert Amber's topology and coordinate files to OpenMM System and Topology stored in XML and PDB files, respectively. 

`setup-settings.sh` includes the list of ligand pairs and their reference alignment atoms, the ATM displacement vector, and the list of residues of the receptor that define the binding site. The script assumes that the `mol2` files of the ligand are stored in the `ligands` subdirectory and the `pdb` file of the receptor is stored in the `receptor` subdirectory. It is assumed that the `pdb` file of the receptor is prepared for Amber (see the [Amber tutorials](https://ambermd.org/tutorials/)). The ligands are assumed to have been docked into the binding site.

The setup creates simulation folders in the `complexes` subdirectory for each RBFE calculation. For example, `eralpha-2d-2e` corresponds to the binding free energy calculation of ligand 2e vs 2d.

After the setup script completes, go to the `complexes` directory to minimize and equilibrate the systems:
```
cd  $HOME/AToM-OpenMM/examples/RBFE/eralpha/complexes
bash ./prep.sh
```
This step prepares the systems in the alchemical intermediate state at λ=1/2. The resulting structures are the input of the alchemical replica exchange simulations.

### Alchemical Replica Exchange

Run replica exchange in each of the simulation folders. For example:
```
cd $HOME/AToM-OpenMM/examples/RBFE/eralpha/complexes
for i in eralpha-* ; do ( cd $i ; bash ./run.sh ) ; done
```
The `run.sh` shell scripts are formatted for optionally running them on a `slurm` queuing system. Edit `run_template.sh` in `$HOME/AToM-OpenMM/examples/RBFE/eralpha/scripts` to adapt them to your cluster.

Each replica exchange calculation is set to run for 2 hours on 1 GPU. Much longer running times (24 hours or more) are needed for this system to approach convergence depending on the speed of the GPU. More GPUs can be deployed by editing the `run_template.sh` files in `$HOME/AToM-OpenMM/examples/RBFE/eralpha/scripts`.

### Free Energy Analysis

The relative binding free energy (ΔΔGb = ΔGb(B) - ΔGb(A)) between each ligand pair A and B are collected by the `free_energies.sh` script in the `complexes` directory:
```
cd $HOME/examples/RBFE/eralpha/complexes
bash ./free_energies.sh
```

The output would look like:
```
2d-2e DDGb=   1.27 +-   0.37  range: 20 35
2d-3a DDGb=  -1.99 +-   0.52  range: 20 35
...
```
etc. where `DDGb` is the value of the relative binding free energy of the second ligand relative to the first. In this example the first 20 samples are discarded and each simulation collected 35 samples per replica. Edit the `free_energies.sh` script to change the number of samples discarded. Obviously much more data is needed to reach convergence for this system in practice. This is only an example.
