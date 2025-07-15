Absolute Binding Free Energy Between the TEMOA host and the G1 guest
--------------------------------------------------------------------

In this tutorial we will calculate the binding free energy of the TEMOA-G1 complex from the [SAMPL8 GDCC](https://github.com/samplchallenges/SAMPL8/tree/master/host_guest/GDCC) challenge using Alchemical Hamiltonian Replica Exchange with the [Alchemical Transfer Method (ATM)](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00266) using [AToM-OpenMM software](https://github.com/Gallicchio-Lab/AToM-OpenMM) and the [ATMMetaForce OpenMM plugin](https://github.com/Gallicchio-Lab/openmm-atmmetaforce-plugin). See [README](https://github.com/Gallicchio-Lab/AToM-OpenMM/blob/master/examples/ABFE/temoa-g1/README.md) for additional specific software requirements.

See [Azimi, Wu, Khuttan, Kurtzman, Deng and Gallicchio.Application of the Alchemical Transfer and Potential of Mean Force Methods to the SAMPL8 Host-Guest Blinded Challenge](https://doi.org/10.1007/s10822-021-00437-y) for further information about the alchemical theory, the ATM method, and this specific application. 

We assume in this tutorial that the OpenMM and AToM-OpenMM packages and their dependencies are available in a `conda` environment (see the [README](../../../README.md), and that the `examples` folder is available under `$HOME/AToM-OpenMM/examples`. Adjust this pathname as needed.

### System preparation

The starting point are the topology and coordinate files of the TEMOA-G1 complex in a water solvent box in the Amber files `temoa-g1.prmtop` and `temoa-g1.inpcrd` provided in this folder. How to prepare systems in Amber format is beyond the scope of this tutorial. We used the `Antechamber` and `tleap` programs of the [`AmberTools` suite version 19](https://ambermd.org/) using the GAFF force field and the TIP3P water model. The Amber topology and coordinate files are converted to OpenMM System XML file and a Topology PDB file.

Prepare, minimize, thermalize, relax, and equilibrate the complex:
```
cd $HOME/AToM-OpenMM/examples/ABFE/temoa-g1
make_atm_system_from_amber --AmberPrmtopinFile temoa-g1.prmtop --AmberInpcrdinFile temoa-g1.inpcrd --systemXMLoutFile temoa-g1_sys.xml --systemPDBoutFile temoa-g1.pdb
abfe_structprep temoa-g1_asyncre.cntl
```

The resulting structure, stored in the `temoa-g1_0.xml` file is the input of the replica exchange production calculation. The PDB version of each structure is also available.

### Replica Exchange

Replica exchange runs on devices listed in a `nodefile`. The provided `nodefile` assumes one GPU on the system. It looks like:
```
localhost,0:0,1,CUDA,,/tmp
```
The critical bit is the `0:0` item in the format `<Platform id>:<device id>`. The other items are in the nodefile specification are for future use and are ignored. You can add more GPUs if you have them. For example, create this nodefile to use two GPUs with device ids 0, and 1:
```
localhost,0:0,1,CUDA,,/tmp
localhost,0:1,1,CUDA,,/tmp
```

To use the OpenCL platform, replace `CUDA` with `OpenCL` in the nodefile, in which case the first number in the device specification should be set to the index of the OpenCL platform corresponding to the GPUs on the system (probably 0). With CUDA, the platform id specification is ignored.

OpenMM's CPU and Hip platforms are also supported. For example,
```
localhost,0:0,1,HIP,,/tmp
localhost-cpu,0:0,4,CPU,,/tmp
```
would use an AMD GPU together with 4 CPU cores on the same system. In most cases, the GPU will be many times faster than the CPU. Nevertheless, the CPU will contribute somewhat to the overall throughput.

Now run replica exchange
```
abfe_production temoa-g1_asyncre.cntl
```

You should see the contents of the control file echo-ed back and messages indicating that replicas are dispatched for execution to the GPUs/CPUs and that replicas change alchemical states by exchanging them with other replicas. 

The job is set to run for 4 hours. The amount of samples collected during this time will depend on the speed of your hardware. To reproduce the equilibrated free energy values in the [paper](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00266) a run of 12 to 24 hours on one GPU would be probably required. The trajectory data for each replica is stored in the `r0`, `r1`, etc. subfolders. These folders contain a `.out` file with one line per sample with perturbation energy and other information, a `.xtc` trajectory file and a checkpoint file to restart the replica exchange simulation. 

The `.out` file of each replica is in the following sample format
```
10 300.000000 -1.000000 0.500000 0.500000 0.100000 30.000000 0.000000 -27850.664985 -10.090523 0.000000
10 300.000000  1.000000 0.500000 0.500000 0.100000 30.000000 0.000000 -27867.121487   4.324510 0.000000
9  300.000000  1.000000 0.400000 0.500000 0.100000 40.000000 0.000000 -28118.427490  24.677342 0.000000
9  300.000000  1.000000 0.400000 0.500000 0.100000 40.000000 0.000000 -28019.335148  27.097275 0.000000
```
where each line is a sample (saved every 5ps in this case), the first column is the alchemical state id (λ=0 is state 0 and λ=1/2 is state 10 in this case), the second column is the set temperature, the third column is the direction state parameter (1 from bound to unbound, and -1 from unbound to bound), the fourth to the to eighth columns hold the softplus alchemical parameters λ1, λ2, α, u0, and w0 (see the [paper](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00266) for details), the 9th column holds the potential energy and the last column is the perturbation energy which is used to compute the free energy change (see below). The 10th column is reserved for future use.

The trajectories can be viewed with VMD. For example this will load the trajectory for replica 3 in VMD:
```
cd $HOME/examples/ABFE/temoa-g1
vmd -f temoa-g1_0.pdb r3/temoa-g1.xtc
```
The output files and the trajectory files can be viewed while replica exchange is running.

Hit `ctrl-C` in the window that runs replica exchange to kill the calculation prematurely. It might take a few seconds for the job to clean up and terminate. The replica exchange job can be restarted by re-issuing the same command as above
```
abfe_production temoa-g1_asyncre.cntl
```
it will restart from the last saved checkpoint. In this example checkpoint files are saved every 10 minutes.

#### Free Energy Analysis

Free energy analysis is performed with [UWHAM](https://cran.r-project.org/web/packages/UWHAM/index.html) in the R statistical environment.

Run the analyze script:
```
./analyze.sh 20
```
It should print something similar to this:
```
DGb = -7.30192 +- 0.2542844 range: 20 38
```
The first number is the estimate of the excess binding free energy in kcal/mol. In this example 20 is the number of initial samples to discard for equilibration. 38 in this case is the number of samples per replica that has been collected. 38 samples with 20 discarded, corresponds in this case to only ~100 ps of data, way too short to obtain a converged value. This is just an example.

The standard binding free energy is obtained by adding the ideal component to the excess component. Again see the [paper](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00266) for details.
