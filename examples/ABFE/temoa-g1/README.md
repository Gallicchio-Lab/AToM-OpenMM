Absolute Binding Free Energy Between The TEMOA host and the G1 guest
--------------------------------------------------------------------

In this tutorial we will calculate the binding free energy of the TEMOA-G1 complex from the [SAMPL8 GDCC](https://github.com/samplchallenges/SAMPL8/tree/master/host_guest/GDCC) challenge using Alchemical Hamiltonian Replica Exchange with the [Alchemical Transfer Method (ATM)](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00266) using [ASyncRE-OpenMM](https://github.com/Gallicchio-Lab/async_re-openmm) and the [ATMMetaForce OpenMM plugin](https://github.com/Gallicchio-Lab/openmm-atmmetaforce-plugin). See [README](https://github.com/Gallicchio-Lab/async_re-openmm/blob/master/examples/ABFE/temoa-g1/README.md) for additional specific software requirements.

See [Azimi, Wu, Khuttan, Kurtzman, Deng and Gallicchio.Application of the Alchemical Transfer and Potential of Mean Force Methods to the SAMPL8 Host-Guest Blinded Challenge](https://arxiv.org/abs/2107.05155) for further information about the alchemical theory, the ATM method, and the chemical systems. 

### System preparation

The starting point are the topology and coordinate files of the TEMOA-G1 complex in a water solvent box in the Amber files `temoa-g1.prmtop` and `temoa-g1.inpcrd` provided in this folder. How to prepare systems in Amber format is beyond the scope of this tutorial. We used the `Antechamber` and `tleap` programs of the [`AmberTools` suite version 19](https://ambermd.org/) using the GAFF force field and the TIP3P water model.

Mininize, thermalize, relax, and equilibrate the complex:
```
python mintherm.py && python npt.py && python equil.py
```
`mintherm` and `npt` equilibrate the solvent keeping the complex restrained. `equil` equilibrates the whole system keeping only the lower cup of the host loosely restrained as in the original work. Each step creates an OpenMM checkpoint file in XML format to start the subsequent step. Each step also generates a PDB file for visualization.

The next step is specific to the ATM method. ATM computes the free energy in two legs that connect the bound and unbound state to the so-called alchemical intermediate that corresponds to an unphysical state which is half bound and half unbound. Check out the [paper](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00266) for more information. The following step slowly ramps up the λ alchemical parameter from zero (corresponding to the bound state) to 1/2 (corresponding to the alchemical intermediate). 
```
python mdlambda.py
```
The resulting structure, stored in the `temoa-g1_0.xml` file is the input of the first leg of the replica exchange calculation. The input of the second leg is stored in the `temoa-g1_0_displaced.xml` file in which the ligand is displaced in the bulk. The PDB version of each structure is also available.

### Replica Exchange

#### Leg 1 - from the bound state to the alchemical intermediate

Copy the input files into the simulation directory
```
cp temoa-g1.prmtop temoa-g1.inpcrd temoa-g1_0.xml asyncre-leg1/
```
Copy also the `nodefile` from the scripts directory
```
cp ../../scripts/nodefile asyncre-leg1/
```
This `nodefile` assumes one GPU on the system (on OpenCL platform 0 with device id 0). It looks like:
```
localhost,0:0,1,OpenCL,,/tmp
```
The critical bit is the `0:0` item in the format `<OpenCL platform id>:<device id>'. The other items are in the nodefile specification are for future use and are ignored. You can add more GPUs if you have them. For example, create this nodefile to use four GPUs with device ids 0, 1, 2, and 3:
```
localhost,0:0,1,OpenCL,,/tmp
localhost,0:1,1,OpenCL,,/tmp
localhost,0:2,1,OpenCL,,/tmp
localhost,0:3,1,OpenCL,,/tmp
```

