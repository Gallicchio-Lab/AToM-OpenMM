RElative Binding Free Energy Between the G1 and G4 Guests for the TEMOA Host
---------------------------------------------------------------------------

In this tutorial we will calculate the difference in the binding free energy between the TEMOA-G4 and TEMOA-G1 complex from the [SAMPL8 GDCC](https://github.com/samplchallenges/SAMPL8/tree/master/host_guest/GDCC) challenge using Alchemical Hamiltonian Replica Exchange with the [Alchemical Transfer Method (ATM)](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00266) using [ASyncRE-OpenMM](https://github.com/Gallicchio-Lab/async_re-openmm) and the [ATMMetaForce OpenMM plugin](https://github.com/Gallicchio-Lab/openmm-atmmetaforce-plugin). See [README](https://github.com/Gallicchio-Lab/async_re-openmm/blob/master/examples/ABFE/temoa-g1/README.md) for additional specific software requirements.

See [Solmaz Azimi, Sheenam Khuttan, Joe Z. Wu, Rajat K. Pal, and Emilio  Gallicchio. Relative Binding Free Energy Calculations for Ligands with Diverse Scaffolds with the Alchemical Transfer Method. ArXiv Preprint 2107.05153 (2021).](https://arxiv.org/abs/2107.05153) for further information about the alchemical theory, the RBFE protocol of the ATM method, and the chemical systems.

### System preparation

The starting point are the topology and coordinate files of a simulation box with the TEMOA host and the G1 and G4 guests in the Amber files `temoa-g1-g4.prmtop` and `temoa-g1-g4.inpcrd` provided in this folder. How to prepare systems in Amber format is beyond the scope of this tutorial. We used the `Antechamber` and `tleap` programs of the [`AmberTools` suite version 19](https://ambermd.org/) using the GAFF force field and the TIP3P water model. G1 was placed in the binding site of the host and G4 at some distance away in the bulk. See the [paper](https://arxiv.org/abs/2107.05153) for more information.

We assume in this tutorial that the examples directory of this repository has been copied under `$HOME/examples` and that the ASyncRE software is available under `$HOME/software/async_re-openmm`. Adjust the pathnames as needed. We are also assuming that OpenMM is launched by running the provided `runopenmm` script. See [examples/README](../../README.md).


Mininize, thermalize, relax, and equilibrate the complex:
```
cd $HOME/RBFE/temoa-g1-g4
../scripts/runopenmm mintherm.py && ../scripts/runopenmm  npt.py && ../scripts/runopenmm equil.py
```
`mintherm` and `npt` equilibrate the solvent keeping the complex restrained. `equil` equilibrates the whole system keeping only the lower cup of the host loosely restrained as in the original work. Each step creates an OpenMM checkpoint file in XML format to start the subsequent step. Each step also generates a PDB file for visualization. All of these steps are performed at the alchemical intermediate state at λ=1/2.

### Replica Exchange Production

```
Copy the `nodefile` from the scripts directory
```
cp ../scripts/nodefile .
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

Now run replica exchange
```
../scripts/runopenmm $HOME/software/async_re-openmm/rbfe_explicit.py temoa-g1-g4_asyncre.cntl
```

You should see the contents of the control file echo-ed back and messages indicating that replica are dispatched to the GPU and that replicas change alchemical states by exchanging them with other replicas. The job is set to run for 4 hours.

#### Free Energy Analysis

Start by copying the necessary scripts to the simulation folder:
```
cd $HOME/examples/RBFE/temoa-g1-g4
cp ../scripts/{analyze.sh,uwham_analysis.R} .
```

Then run the analyze script:
```
./analyze.sh 20
```
It should print something similar to this:
```
DDGb = -3.50378 +- 0.3484641 range: 20 39
```
The value shown is the relative binding free energy (ΔΔGb = ΔGb(TEMOA-G4) - ΔG(TEMOA-G1)) between TEMOA-G4 and TEMOA-G1. In this example 20 is the number of initial samples to discard for equilibration. 39 in this case is the number of samples per replica that has been collected. The simulations here are just short examples. See the [paper](https://arxiv.org/abs/2107.05153) for more details.
