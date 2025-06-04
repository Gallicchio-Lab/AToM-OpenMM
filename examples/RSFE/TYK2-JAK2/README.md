Relative Binding Selectivity Free Energy of a pair of ligands to TYK2 and JAK2 
------------------------------------------------------------------------------

In this example, a swapping calculation will be performed to calculate the relative binding selectivity free energy between a pair of ligands to the isoforms TYK2 and JAK2 of the JAK tyrosine kinase family (Liang et al. (2013): Lead identification of novel and selective TYK2 inhibitors, https://pubmed.ncbi.nlm.nih.gov/23867602/). This example uses an automated workflow to prepare and run the calculation. 


### System preparation

Assuming AToM-OpenMM is available under `$HOME/AToM-OpenMM`. Input the pathname that goes to the AToM-OpenMM installation folder in  `scripts/setup-settings.sh`. The setup script `setup-atm.sh` reads parameters for the system from the  `scripts/setup-settings.sh` file. 

```
cd $HOME/AToM-OpenMM/examples/RSFE/TYK2-JAK2
bash ./scripts/setup-atm.sh
```

The setup scripts uses OpenMM and `openmmforcefields`  to prepare the ligands and receptors. For the ligands, the OpenFF-2.0.0 force field is used and for the receptors, the protein Amber14 forcefield. The system is solvated in a box. `setup-settings.sh` includes the ligand pairs and their reference alignment atoms. For swapping, alignment atoms are not required, however they are necessary inputs for the workflow. If RBFE calculations are sought to compliment swapping calculations, the alignment atoms selected in this example are sound. Additional ligands can be included and their details must be inputted in this script to automatically set up additional systems.

The displacement vector is necessary for swapping. A sufficiently large displacement must be selected to ensure at least three layers of water molecules between the two proteins. The setup script assumes that the `sdf` files of the ligand are stored in the `ligands` subdirectory and the `pdb` file of the receptor is stored in the `receptor` subdirectory. It is assumed that the pdb file of the receptor is fully prepared with hydrogen atoms etc. The ligands are assumed to have been docked into the binding site.

The setup creates simulation folders in the `complexes` subdirectory for each swapping calculation. For example, `TYK2-JAK2-c41-c46` corresponds to the swapping of c41 from TYK2 to JAK2 and the subsequent swap of c46 from JAK2 to TYK2. 

### Minimization and Alchemical Replica Exchange

After setting up the simulation folders, enter the `complexes` subdirectory to run minimization and equilibration for each system. These processes prepare the systems at the alchemical intermediate state at λ=1/2 and the resulting configurations are inputs for alchemical replica exchange simulations. Then alchemical replica exchange can be conducted for each system. 

```
cd $HOME/AToM-OpenMM/examples/RSFE/TYK2-JAK2/complexes
for i in TYK2-JAK2-* ; do ( cd $i ; bash ./run.sh ) ; done
```

The `run.sh` shell scripts are formatted for a slurm queuing system. Edit run_template.sh in `$HOME/AToM-OpenMM/examples/RSFE/TYK2-JAK2/scripts` to adapt them to the cluster of choice. 

Each simulation is set to run for 24 hours on 1 GPU. Swapping calculations often require longer run times to reach convergence. More GPUs can be deployed by editing the `run_template.sh` file.


### Free Energy Analysis

The relative binding selectivity free energy between the ligand pair A and B are collected by the `free_energies.sh` script in the complexes directory:
```
cd $HOME/AToM-OpenMM/examples/RSFE/TYK2-JAK2/complexes
bash ./free_energies.sh
```

The `free_energies.sh` script produces the swapping free energy estimate and the statistical uncertainty associated with it. In addition, it will indicate if there are large gaps in the distributions of the binding free energies, which would indicate the calculation has not converged. A plot of the distributions, as well as the free energy profile, are also produced for further analysis.

Edit the `free_energies.sh` script to change the number of samples discarded. It is customary to discard the first half of total samples.

If it throws an error named `./analyze.sh: line 34: R: command not found`, please install `r-base` and `UWHAM R package` as follow [installation instructions](https://github.com/Gallicchio-Lab/AToM-OpenMM#installation--usage).

