RElative Binding Free Energy Between the G1 and G4 Guests for the TEMOA Host
---------------------------------------------------------------------------

In this tutorial we will calculate the difference in the binding free energy between the TEMOA-G4 and TEMOA-G1 complex from the [SAMPL8 GDCC](https://github.com/samplchallenges/SAMPL8/tree/master/host_guest/GDCC) challenge using Alchemical Hamiltonian Replica Exchange with the [Alchemical Transfer Method (ATM)](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00266) using [ASyncRE-OpenMM](https://github.com/Gallicchio-Lab/async_re-openmm) and the [ATMMetaForce OpenMM plugin](https://github.com/Gallicchio-Lab/openmm-atmmetaforce-plugin). See [README](https://github.com/Gallicchio-Lab/async_re-openmm/blob/master/examples/ABFE/temoa-g1/README.md) for additional specific software requirements.

See [Solmaz Azimi, Sheenam Khuttan, Joe Z. Wu, Rajat K. Pal, and Emilio  Gallicchio. Relative Binding Free Energy Calculations for Ligands with Diverse Scaffolds with the Alchemical Transfer Method. ArXiv Preprint 2107.05153 (2021).](https://arxiv.org/abs/2107.05153) for further information about the alchemical theory, the RBFE protocol of the ATM method, and the chemical systems.

### System preparation

The starting point are the topology and coordinate files of a simulation box with the TEMOA host and the G1 and G4 guests in the Amber files `temoa-g1-g4.prmtop` and `temoa-g1-g4.inpcrd` provided in this folder. How to prepare systems in Amber format is beyond the scope of this tutorial. We used the `Antechamber` and `tleap` programs of the [`AmberTools` suite version 19](https://ambermd.org/) using the GAFF force field and the TIP3P water model. G1 was placed in the binding site of the host and G4 at some distance away in the bulk. See the [paper](https://arxiv.org/abs/2107.05153) for more information.

We assume in this tutorial that the examples folder of this repo has been copied under `$HOME/examples`. Adjust the pathnames as needed.

Mininize, thermalize, relax, and equilibrate the complex:
```
cd $HOME/RBFE/temoa-g1-g4
python mintherm.py && python npt.py && python equil.py
```
`mintherm` and `npt` equilibrate the solvent keeping the complex restrained. `equil` equilibrates the whole system keeping only the lower cup of the host loosely restrained as in the original work. Each step creates an OpenMM checkpoint file in XML format to start the subsequent step. Each step also generates a PDB file for visualization. All of these steps are performed at the alchemical intermediate state at λ=1/2.

### Replica Exchange

#### Leg 1 - from G1 bound to G4 unbound to the alchemical intermediate

Copy the input files into the simulation directory
```
cp temoa-g1-g4.prmtop temoa-g1-g4.inpcrd temoa-g1-g4_0.xml asyncre-leg1/
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
python rbfe_explicit.py temoa-g1-g4_asyncre.cntl
```

You should see the contents of the control file echo-ed back and messages indicating that replica are dispatched to the GPU and that replicas change alchemical states by exchanging them with other replicas. The job is set to run for two hours.

#### Leg 2 - from the unbound state to the alchemical intermediate

The leg 2 calculation can run in parallel to the leg 1 calculation if you have multiple computing nodes or GPU devices available. Copy the input files to the replica exchange folder for leg 2:
```
cd $HOME/examples/RBFE/temoa-g1-g4
cp temoa-g1-g4.prmtop temoa-g1-g4.inpcrd asyncre-leg2/
cp temoa-g1-g4_0_displaced.xml asyncre-leg2/temoa-g1-g4_0.xml
cp ../../scripts/nodefile asyncre-leg2/
```
Notice that this time we copied the structure with the ligand displaced. Then run replica exchange as before

Now go the replica exchange folder for leg 2 and run replica exchange
```
cd asyncre-leg2/
python rbfe_explicit.py temoa-g1-g4_asyncre.cntl
```


