RElative Binding Free Energy Between the G1 and G4 Guests for the TEMOA Host
---------------------------------------------------------------------------

In this tutorial we will calculate the difference in the binding free energy between the TEMOA-G4 and TEMOA-G1 complex from the [SAMPL8 GDCC](https://github.com/samplchallenges/SAMPL8/tree/master/host_guest/GDCC) challenge using Alchemical Hamiltonian Replica Exchange with the [Alchemical Transfer Method (ATM)](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00266) using [AToM-OpenMM](https://github.com/Gallicchio-Lab/AToM-OpenMM) and the [ATM MetaForce OpenMM plugin](https://github.com/Gallicchio-Lab/openmm-atmmetaforce-plugin). See the [README](https://github.com/Gallicchio-Lab/AToM-OpenMM/blob/master/examples/ABFE/temoa-g1/README.md) for additional specific software requirements.

See [Solmaz Azimi, Sheenam Khuttan, Joe Z. Wu, Rajat K. Pal, and Emilio  Gallicchio. Relative Binding Free Energy Calculations for Ligands with Diverse Scaffolds with the Alchemical Transfer Method](https://pubs.acs.org/doi/10.1021/acs.jcim.1c01129) for further information about the alchemical theory, the RBFE protocol of the ATM method, and the chemical systems.

It is highly recommended to go through the [TEMOA G1 ABFE tutorial](https://github.com/Gallicchio-Lab/AToM-OpenMM/tree/master/examples/ABFE/temoa-g1) before attempting this one. The following assumes familiarity with the terms and procedures introduced in the [TEMOA G1 ABFE tutorial](https://github.com/Gallicchio-Lab/AToM-OpenMM/tree/master/examples/ABFE/temoa-g1).

We assume in this tutorial that the OpenMM and AToM-OpenMM packages and their dependencies are available in a `conda` environment (see the [README](../../../README.md), and that the `examples` folder is available under `$HOME/AToM-OpenMM/examples`. Adjust this pathname as needed.

### System preparation

The starting point are the topology and coordinate files of a simulation box with the TEMOA host and the G1 and G4 guests in the Amber files `temoa-g1-g4.prmtop` and `temoa-g1-g4.inpcrd` provided in this folder. How to prepare systems in Amber format is beyond the scope of this tutorial. We used the `Antechamber` and `tleap` programs of the [`AmberTools` suite version 19](https://ambermd.org/) using the GAFF force field and the TIP3P water model. G1 was placed in the binding site of the host and G4 at some distance away in the bulk. See the [paper](https://pubs.acs.org/doi/10.1021/acs.jcim.1c01129) for more information.

Prepare, minimize, thermalize, relax, and equilibrate the complex:
```
cd $HOME/AToM-OpenMM/examples/RBFE/temoa-g1-g4
make_atm_system_from_amber --AmberPrmtopinFile temoa-g1-g4.prmtop --AmberInpcrdinFile temoa-g1-g4.inpcrd --systemXMLoutFile temoa-g1-g4_sys.xml --systemPDBoutFile temoa-g1-g4.pdb
rbfe_structprep temoa-g1-g4_asyncre.cntl
```
The lower cup of the host loosely restrained as in the original work. Each step creates an OpenMM checkpoint file in XML format to start the subsequent step. Each step also generates a PDB file for visualization. The result is an equilibrated system at the alchemical intermediate state at λ=1/2.

### Replica Exchange Production

See the [TEMOA-G1 ABFE tutorial](https://github.com/Gallicchio-Lab/AToM-OpenMM/tree/master/examples/ABFE/temoa-g1) about the nodefile and how to customize it in different ways to match the hardware on your machine.

```
rbfe_production temoa-g1-g4_asyncre.cntl
```

You should see the contents of the control file echo-ed back and messages indicating that replica are dispatched to the GPU and that replicas change alchemical states by exchanging them with other replicas. The job is set to run for 4 hours.

#### Free Energy Analysis

Run the analyze script:
```
./analyze.sh 20
```
It should print something similar to this:
```
DDGb = -3.50378 +- 0.3484641 range: 20 39
```
The value shown is the relative binding free energy (ΔΔGb = ΔGb(TEMOA-G4) - ΔGb(TEMOA-G1)) between TEMOA-G4 and TEMOA-G1. In this example 20 is the number of initial samples to discard for equilibration. 39 in this case is the number of samples per replica that has been collected. The simulations here are just short examples. See the [paper](https://pubs.acs.org/doi/10.1021/acs.jcim.1c01129) for more details.
