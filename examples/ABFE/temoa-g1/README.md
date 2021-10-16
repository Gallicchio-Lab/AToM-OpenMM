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
`mintherm` and `npt` equilibrate the solvent keeping the complex restrained. `equil` equilibrates the whole system keeping only the lower cup of the host loosely restrained.
