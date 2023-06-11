Absolute Binding Free Energies of a Set of Complexes of the FKBP Protein Receptor with Ligand Fragments 
-------------------------------------------------------------------------------------------------------

In this tutorial we will calculate the absolute binding free energies of the complexes of FKBP with 7 ligand fragments from the work of [A Pan, H Xu, T Palpant and DE Shaw](http://dx.doi.org/10.1021/acs.jctc.7b00172). This tutorial uses an automated workflow to prepare and run the 7 calculations.

It is highly recommended to go through the [TEMOA G1 ABFE tutorial](https://github.com/Gallicchio-Lab/AToM-OpenMM/tree/master/examples/ABFE/temoa-g1) before attempting this one. The following assumes familiarity with the terms and procedures introduced in the [TEMOA G1 ABFE tutorial](https://github.com/Gallicchio-Lab/AToM-OpenMM/tree/master/examples/ABFE/temoa-g1).

### System preparation

We assume in this tutorial that this repository has been cloned under `$HOME`. Adjust the pathname pointing to the AToM-OpenMM installation folder in `setup-settings.sh` as needed. See the [TEMOA G1 ABFE tutorial](https://github.com/Gallicchio-Lab/AToM-OpenMM/tree/master/examples/ABFE/temoa-g1) for help setting up OpenMM and the AToM-OpenMM software.

This tutorial assumes that the [Ambertools](http://ambermd.org) executables are in the search path. They can be conveniently installed under the same conda environment:
```
conda install -c conda-forge ambertools
```

Setup the simulation input files. The automated setup script below reads the parameters from the `setup-settings.sh` file.
```
cd $HOME/AToM-OpenMM/examples/ABFE/fkbp/
bash ./scripts/setup-atm.sh
```
The automated script uses OpenMM's Modeller class, openmmforcefields, and the OpenFF utilities to prepare the ligands and the receptor and to solvated them in a solution box. `setup-settings.sh` includes the list of ligands, the ATM displacement vector, and the list of residues of the receptor that define the binding site. The script assumes that the `sdf` files of the ligand are stored in the `ligands` subdirectory and the `pdb` file of the receptor is stored in the `receptor` subdirectory. It is assumed that the `pdb` file of the receptor is prepared for Amber (see the [Amber tutorials](https://ambermd.org/tutorials/)). The ligands are assumed to have been docked into the binding site.

The setup creates simulation folders in the `complexes` subdirectory for each ABFE calculation. For example, `fkbp-dss` corresponds to the binding free energy calculation for the complex of FKBP with ligand dss.

After the setup script completes, go to the `complexes` directory to minimize and equilibrate the systems:
```
cd  $HOME/examples/ABFE/fkbp/complexes
bash ./prep.sh
```
This step prepares the systems in the alchemical intermediate state at λ=1/2. The resulting structures are the input of the alchemical replica exchange simulations.

### Alchemical Replica Exchange

Run replica exchange in each of the simulation folders. For example:
```
cd $HOME/AToM-OPenMM/examples/ABFE/fkbp/complexes
for i in fkbp-* ; do ( cd $i ; bash ./run.sh ) ; done
```
The `run.sh` shell scripts are formatted for optionally running them on a `slurm` queuing system. Edit `run_template.sh` in `$HOME/examples/ABFE/fkbp/scripts` to adapt them to your cluster.

Each replica exchange calculation is set to run for 2 hours on 1 GPU. Much longer running times (24 hours or more) are needed for this system to approach convergence depending on the speed of the GPU. More GPUs can be deployed by editing the `run_template.sh` file in `$HOME/examples/ABFE/fkbp/scripts`.

### Free Energy Analysis

The binding free energies (ΔGb) of each complex are computed by the `free_energies.sh` script in the `complexes` directory:
```
cd $HOME/examples/ABFE/fkbp/complexes
bash ./free_energies.sh
```

The output would look like:
```
fkbp-thi DGb=  -3.18 +-   0.37  range: 20 35
fkbo-dss DGb=  -2.99 +-   0.52  range: 20 35
...
```
etc. where `DGb` is the value of the binding free energy. In this example the first 20 samples are discarded and each simulation collected 35 samples per replica. Edit the `free_energies.sh` script to change the number of samples discarded. Obviously much more data is needed to reach convergence for this system in practice. This is only an example.
