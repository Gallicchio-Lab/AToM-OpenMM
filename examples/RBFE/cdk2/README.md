Relative Binding Free Energies of a set of ligands of the CDK2 Kinase
---------------------------------------------------------------------


In this tutorial we will calculate the relative binding free energies between congeneric ligands binding to the CDK2 receptor with six ligands from the paper: [Wang et al. (2013) Modeling Local Structural Rearrangements Using FEP/REST: Application to Relative Binding Affinity Predictions of CDK2 Inhibitors](https://pubs.acs.org/doi/10.1021/ct300911a). This tutorial uses an automated workflow to prepare and run the eight relative binding free energy calculations between all pairs of complexes.

We assume that the OpenMM and AToM-OpenMM packages and their dependencies are available in a `conda` environment (see the [README](../../../README.md), and that the `examples` folder is available under `$HOME/AToM-OpenMM/examples`. Adjust this pathname as needed.

### System preparation

Setup the simulation input files. The automated setup script below reads the parameters from the `setup-settings.sh` file.

```
cd $HOME/AToM-OpenMM/examples/RBFE/cdk2
bash ./scripts/setup-atm.sh
```
The automated script uses OpenMM and `openmmforcefields` to prepare the ligands and the receptor using the Sage OpenFF force field and protein Amber force fields and to solvate them in a solution box.

`setup-settings.sh` includes the list of ligand pairs and their reference alignment atoms chosen along the xyz coordinates of the ligand-pair, the ATM displacement vector (to displace the ligand at 
the nearest distance away from the receptor), and the list of residues of the receptor that define the binding site. For more informaion, refer to [ Azimi, Khuttan, Wu, Pal, Gallicchio. Relative Binding Free Energy Calculations for Ligands with Diverse Scaffolds with the Alchemical Transfer Method.](https://pubs.acs.org/doi/full/10.1021/acs.jcim.1c01129)

The script assumes that the `sdf` files of the ligand are stored in the `ligands` subdirectory and the `pdb` file of the receptor is stored in the `receptor` subdirectory. It is assumed that the pdb file of the receptor is fully prepared with hydrogen atoms etc. The ligands are assumed to have been docked into the binding site.

The setup creates simulation folders in the complexes subdirectory for each RBFE calculation. For example, `cdk2-H1Q-OIU` corresponds to the binding free energy calculation of ligand `OIU` vs `H1Q`

### Minimization and Alchemical Replica Exchange

After the setup script completes, go to the `complexes` directory to run minimization and equilibration followed by Alchemical replica exchange for each of the complexes. Minimization and equilibration steps prepare the systems in the alchemical intermediate state at λ=1/2. The resulting structures are the input of the alchemical replica exchange simulations.

For example:
```
cd $HOME/examples/RBFE/cdk2/complexes
for i in cdk2-* ; do ( cd $i ; bash ./run.sh ) ; done
```
The `run.sh` shell scripts are formatted for optionally running them on a slurm queuing system, such as the one on the Expanse supercomputer. Edit run_template.sh in `$HOME/AToM-OpenMM/examples/RBFE/cdk2/scripts` to adapt them to your cluster.

Each replica exchange calculation is set to run for 8 hours on 1 GPU. Much longer running times (24 hours or more) are needed for this system to approach convergence depending on the speed of the GPU. More GPUs can be deployed by editing the `run_template.sh` file in `$HOME/examples/RBFE/cdk2/scripts`

### Free Energy Analysis
The relative binding free energy (ΔΔGb = ΔGb(B) - ΔGb(A)) between each ligand pair A and B are collected by the `free_energies.sh` script in the complexes directory:

```
cd $HOME/AToM-OpenMM/examples/RBFE/cdk2/complexes
bash ./free_energies.sh
```

The free energy values produced here may or may not match with the literature set mentioned in the [reference paper](https://pubs.acs.org/doi/10.1021/ct300911a). Some ligands, `1H1R` and `1H1S` in particular, need enhanced sampling to reach convergence. The mentioned ligands are trapped in their starting conformations in the binding site failing to adopt alternative conformations because of certain torsional barriers. 

Edit the `free_energies.sh` script to change the number of samples discarded. Usually the first half of total samples are disarded. 

If it throws an error named `./analyze.sh: line 34: R: command not found`, please install `r-base` and `UWHAM R package` as follow [installation instructions](https://github.com/Gallicchio-Lab/AToM-OpenMM#installation--usage).
