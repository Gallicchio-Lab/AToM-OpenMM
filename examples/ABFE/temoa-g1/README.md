Absolute Binding Free Energy Between The TEMOA host and the G1 guest
--------------------------------------------------------------------

In this tutorial we will calculate the binding free energy of the TEMOA-G1 complex from the [SAMPL8 GDCC](https://github.com/samplchallenges/SAMPL8/tree/master/host_guest/GDCC) challenge using Alchemical Hamiltonian Replica Exchange with the [Alchemical Transfer Method (ATM)](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00266) using [ASyncRE-OpenMM](https://github.com/Gallicchio-Lab/async_re-openmm) and the [ATMMetaForce OpenMM plugin](https://github.com/Gallicchio-Lab/openmm-atmmetaforce-plugin). See [README](https://github.com/Gallicchio-Lab/async_re-openmm/blob/master/examples/ABFE/temoa-g1/README.md) for additional specific software requirements.

See [Azimi, Wu, Khuttan, Kurtzman, Deng and Gallicchio.Application of the Alchemical Transfer and Potential of Mean Force Methods to the SAMPL8 Host-Guest Blinded Challenge](https://arxiv.org/abs/2107.05155) for further information about the alchemical theory, the ATM method, and the chemical systems. 

### System preparation

The starting point are the topology and coordinate files of the TEMOA-G1 complex in a water solvent box in the Amber files `temoa-g1.prmtop` and `temoa-g1.inpcrd` provided in this folder. How to prepare systems in Amber format is beyond the scope of this tutorial. We used the `Antechamber` and `tleap` programs of the [`AmberTools` suite version 19](https://ambermd.org/) using the GAFF force field and the TIP3P water model.

We assume in this tutorial that this ABFE folder of the examples directory has been copied under `$HOME/ABFE`. Adjust the pathnames as needed.

Mininize, thermalize, relax, and equilibrate the complex:
```
cd $HOME/ABFE/temoa-g1
python mintherm.py && python npt.py && python equil.py
```
`mintherm` and `npt` equilibrate the solvent keeping the complex restrained. `equil` equilibrates the whole system keeping only the lower cup of the host loosely restrained as in the original work. Each step creates an OpenMM checkpoint file in XML format to start the subsequent step. Each step also generates a PDB file for visualization.

The next step is specific to the ATM method. ATM computes the free energy in two legs that connect the bound and unbound state to the so-called alchemical intermediate that corresponds to an unphysical state which is half bound and half unbound. Check out the [paper](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00266) for more information. The following step slowly ramps up the λ alchemical parameter from zero (corresponding to the bound state) to 1/2 (corresponding to the alchemical intermediate). 
```
python mdlambda.py
```
The resulting structure, stored in the `temoa-g1_0.xml` file is the input of the first leg of the replica exchange calculation. The input of the second leg is stored in the `temoa-g1_0_displaced.xml` file in which the ligand is displaced in the bulk. The PDB version of each structure is also available.

### Replica Exchange

#### Leg 1 - From the bound state to the alchemical intermediate

Copy the input files into the simulation directory
```
cp temoa-g1.prmtop temoa-g1.inpcrd temoa-g1_0.xml asyncre-leg1/
```
Copy also the `nodefile` from the scripts directory
```
cp ../../scripts/nodefile asyncre-leg1/
```
This `nodefile` assumes one GPU on the system (on the OpenCL platform 0 with device id 0). It looks like:
```
localhost,0:0,1,OpenCL,,/tmp
```
The critical bit is the `0:0` item in the format `<OpenCL platform id>:<device id>'. The other items are in the nodefile specification are for future use and are ignored. You can add more GPUs if you have them. For example, create this nodefile to use two GPUs with device ids 0, and 1:
```
localhost,0:0,1,OpenCL,,/tmp
localhost,0:1,1,OpenCL,,/tmp
```

Now go the replica exchange folder for leg 1 and run replica exchange
```
cd asyncre-leg1/
python abfe_explicit.py temoa-g1_asyncre.cntl
```

You should see the contents of the control file echo-ed back and messages indicating that replica are dispatched to the GPU and that replicas change alchemical states by exchanging them with other replicas. 

The job is set to run for two hours. The amount of samples collected during this time will depend on the speed of your GPU. To reproduce the equilibrated free energy values in the [paper](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00266) a run of 12 to 24 hours on one GPU would be probably required. The trajectory data for each replica is stored in the `r0`, `r1`, etc. subfolders. These folders contain a `.out` file with one line per sample with perturbation energy and other information, a `.dcd` trajectory file and a checkpoint file to restart the replica exchange simulation. 

The `.out` file of each replica is in the following sample format
```
10 300.000000 0.500000 0.500000 0.100000 30.000000 0.000000 -27850.664985 10.090523
10 300.000000 0.500000 0.500000 0.100000 30.000000 0.000000 -27867.121487 4.324510
9 300.000000 0.400000 0.500000 0.100000 40.000000 0.000000 -28118.427490 24.677342
9 300.000000 0.400000 0.500000 0.100000 40.000000 0.000000 -28019.335148 27.097275
```
where each line is a sample (saved every 5ps in this case), the first column is the alchemical state id (λ=0 is state 0 and λ=1/2 is state 10 in this case), the second column is the set temperature, the third to the to seventh column hold the softplus alchemical parameters λ1, λ2, α, u0, and w0 (see the [paper](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00266) for details), the 8th column holds the potential energy and the last column is the perturbation energy which is used to compute the free energy change (see below).

The trajectories can be viewed with VMD. For example this will load the trajectory for replica 3 in VMD:
```
cd $HOME/ABFE/temoa-g1
vmd -f temoa-g1_0.pdb asyncre-leg1/r3/temoa-g1.dcd
```
The output files and the trajectory files can be viewed while replica exchange is running.

Hit `ctrl-C` in the window that runs replica exchange to kill the calculation prematurely. It might take a few seconds for the job to clean up and terminate. The replica exchange job can be restarted by re-issuing the same command as above
```
python abfe_explicit.py temoa-g1_asyncre.cntl
```
it will restart from the last saved checkpoint. In this example checkpoint files are saved every 10 minutes.

#### Leg 2 - From the unbound state to the alchemical intermediate

The leg 2 calculation can run in parallel to the leg 1 calculation if you have multiple computing nodes or GPU devices available. Copy the input files to the replica exchange folder for leg 2:
```
cd $HOME/ABFE/temoa-g1
cp temoa-g1.prmtop temoa-g1.inpcrd asyncre-leg2/
cp temoa-g1_0_displaced.xml asyncre-leg2/temoa-g1_0.xml
cp ../../scripts/nodefile asyncre-leg2/
```
Notice that this time we copied the structure with the ligand displaced. Then run replica exchange as before
```
cd asyncre-leg2/
python abfe_explicit.py temoa-g1_asyncre.cntl
```

#### Free Energy Analysis


