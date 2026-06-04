# AToM-OpenMM

The Alchemical Transfer Method for OpenMM (AToM-OpenMM) is an extensible Python package for estimating absolute and relative binding free energies of molecular complexes.

AToM-OpenMM implements the Alchemical Transfer Method (ATM) with asynchronous parallel replica exchange molecular dynamics using [OpenMM](https://github.com/openmm). It can be deployed on workstations or cluster nodes with one or more GPUs.

This version uses [ATMForce](https://github.com/openmm/openmm/pull/4110), available in OpenMM 8.4.0 and later.

## Start Here

- [Installation](installation.md): install AToM-OpenMM and its runtime dependencies.
- [Tutorials](tutorials/index.md): start with Colab notebooks or local repository tutorials.
- [User guide](user-guide/index.md): learn the high-level ABFE and RBFE workflow layout.
- [Theory](theory/index.md): introduces the theory of the Alchemical Transfer approach
- [API reference](api-reference/index.md): find command-line entry points, workflow modules, and utility modules.

## Credits

AToM-OpenMM is developed and maintained by [Emilio Gallicchio's lab](http://www.compmolbiophysbc.org), with support from current and past grants from the National Science Foundation and the National Institutes of Health.

Maintainer and author:

- Emilio Gallicchio

Contributors:

- Elian Tiudic
- Sylvester Sakyi
- Stefan Doerr
- Sheenam Khuttan
- Joe Z Wu
- Solmaz Azimi
- Baofeng Zhang
- Rajat Pal

The asynchronous replica exchange method was first implemented in the [AsyncRE](https://github.com/ComputationalBiophysicsCollaborative/AsyncRE) package for the IMPACT program.

## Citations

Please cite [AToM-OpenMM](http://www.compmolbiophysbc.org/publications) if you use this software in your research.

- [Relative Binding Free Energy Calculations for Ligands with Diverse Scaffolds with the Alchemical Transfer Method](https://pubs.acs.org/doi/10.1021/acs.jcim.1c01129)
- [Relative Binding Free Energy Estimation of Congeneric Ligands and Macromolecular Mutants with the Alchemical Transfer with Coordinate Swapping Method](https://doi.org/10.1021/acs.jcim.5c00207)
- [Alchemical Transfer Approach to Absolute Binding Free Energy Estimation](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00266)
- [Asynchronous Replica Exchange Software for Grid and Heterogeneous Computing](http://www.compmolbiophysbc.org/publications#asyncre_software_2015)

## License

AToM-OpenMM is research software under active development. It is provided without guarantees of correctness.

The software is licensed under the GNU Lesser General Public License. See the repository [LICENSE](https://github.com/Gallicchio-Lab/AToM-OpenMM/blob/master/LICENSE) file for details. The AToM logo is copyright 2023 Solmaz Azimi.
