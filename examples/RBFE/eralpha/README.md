Relative Binding Free Energies of a set of Antagonists of the Estrogen Receptor α Nuclear Receptor
--------------------------------------------------------------------------------------------------

In this tutorial we will calculate the relative binding free energies between four protein-ligand complexes of the ERα nuclear receptor with four antagonists from the paper: [Azimi, Khuttan, Wu, Pal, Gallicchio. Relative Binding Free Energy Calculations for Ligands with Diverse Scaffolds with the Alchemical Transfer Method.](https://arxiv.org/abs/2107.05153) This tutorial uses an automated workflow to prepare and run the six relative binding free energy calculations.

### System preparation

We assume in this tutorial that the examples folder of this repo has been copied under `$HOME/examples`. Adjust the pathnames as needed. Consult this [README](https://github.com/Gallicchio-Lab/async_re-openmm/blob/master/examples/README.md) for software requirements.

Setup the simulation input files. The automated setup script below reads the parameters from the `setup-settings.sh` file. Modify the `work_dir` setting there if the simulation folder is different than the one assumed (`$HOME/examples/RBFE/asyncre`).
```
cd $HOME/examples/RBFE/asyncre/scripts
bash ./setup-atm.sh
```
The automated script runs `antechamber` and `tleap` from the AmberTools package to prepare the ligands and the receptor and to solvated them in a solution box. `setup-settings.sh` also include the list of ligand pairs and their reference alignment atoms, the ATM displacement vector, and the list of residues of the receptor that define the binding site. The script assumes that the `mol2` files of the ligand are stored in the `ligands` subdirectory and the `pdb` file of the receptor is stored in the `receptor` subdirectory. It is assumed that the `pdb` file of the receptor is prepared for Amber (see the [Amber tutorials](https://ambermd.org/tutorials/)). The ligands are assumed placed into the binding site.

The setup creates simulation folders in the `complexes` subdirectory. There are two folders for each ligand pair. For example `eralpha-2d-2e` corresponds to leg 1 for the relative binding free energy calculation of ligand 2e vs 2d. The folder `eralpha-2e-2d` corresponds to leg 2 of the same ligand pair.

After the setup script completes, go to the `complexes` directory to minimize and equilibrate the systems:
```
cd  $HOME/examples/RBFE/asyncre/complexes
bash ./prep.sh
```
This step prepares the systems in the alchemical intermediate state at λ=1/2. The resulting structures are the input of the alchemical replica exchange simulations.

### Alchemical Replica Exchange

Run replica exchange in each of the simulation folders. For example:
```
cd $HOME/examples/RBFE/asyncre/complexes
for i in eralpha-* ; do ( cd $i ; bash ./run.sh ) ; done
```
The `run.sh` shell scripts are formatted for optionally running them on a `slurm` queuing system. Edit `run_template.sh` in `$HOME/examples/RBFE/asyncre/scripts` to adapt them to your cluster.

Each replica exchange calculation is set to run for 2 hours on 1 GPU. Much longer running times (24 hours or more) are needed to approach convergence depending on the speed of the GPU. More GPUs can be deployed by editing the `nodefile` or the `run_template.sh` files in `$HOME/examples/RBFE/asyncre/scripts`.

### Free Energy Analysis

The relative binding free energy (ΔΔGb = ΔGb(B) - ΔG(A)) between each ligand pair A and B is the free energy change in the corresponding leg 1 minus the free energy change in leg 2. Run the automated script in the `complexes` directory:
```
cd $HOME/examples/RBFE/asyncre/complexes
bash ./free_energies.sh
```

The output looks like:
```
2d-2e   DGb = 20.01217 +- 0.1184765 DE = 0.8511286 +- 0.2338686  range: 100 600
2e-2d   DGb = 18.74483 +- 0.1176946 DE = -0.4098593 +- 0.230605  range: 100 600
2d-2e DDGb=   1.27 +-   0.17 DDE =    0.44 +-   0.33 
2d-3a   DGb = 21.76995 +- 0.1409262 DE = -1.070465 +- 0.2560702  range: 100 600
3a-2d   DGb = 23.7594 +- 0.1405318 DE = 1.320733 +- 0.2255252  range: 100 600
2d-3a DDGb=  -1.99 +-   0.20 DDE =    0.25 +-   0.34 
```
etc. Where the first 2 lines report the free energy changes in each leg and the third the relative binding free energy. And so on for each ligand pair.
