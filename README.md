AToM-OpenMM v8.1.2
====================

The Alchemical Transfer Method for OpenMM (AToM-OpenMM) is an extensible Python package for the estimation of absolute and relative binding free energies of molecular complexes. It implements the [Alchemical Transfer Method (ATM)](https://pubs.acs.org/doi/10.1021/acs.jcim.1c01129) with  asynchronous parallel replica exchange molecular dynamics with the [OpenMM](https://github.com/openmm) library. The AToM software can be deployed on workstations or cluster nodes with one or more GPUs.

This version of AToM uses [ATMForce](https://github.com/openmm/openmm/pull/4110) in the 8.1.1 or later versions of [OpenMM](https://github.com/openmm/openmm). 

Credits
-------

This software is developed and maintained by the [Emilio Gallicchio's lab](http://www.compmolbiophysbc.org) with support from current and past grants from the National Science Foundation (ACI 1440665 and CHE 1750511) and the National Institutes of Health (R15 GM151708).

Maintainers:

- Emilio Gallicchio <egallicchio@brooklyn.cuny.edu>

Authors:

- Emilio Gallicchio <egallicchio@brooklyn.cuny.edu>
- Baofeng Zhang <BZhang@brooklyn.cuny.edu>
- Rajat Pal <rajatfor2014@gmail.com>

The asynchronous replica exchange method was first implemented in the [AsyncRE](https://github.com/ComputationalBiophysicsCollaborative/AsyncRE) package for the IMPACT program.

Citations
---------

Please [cite us](http://www.compmolbiophysbc.org/publications) if you use this software in your research:

- [Relative Binding Free Energy Calculations for Ligands with Diverse Scaffolds with the Alchemical Transfer Method](https://pubs.acs.org/doi/10.1021/acs.jcim.1c01129)

- [Alchemical Transfer Approach to Absolute Binding Free Energy Estimation](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00266)

- [Asynchronous Replica Exchange Software for Grid and Heterogeneous Computing](http://www.compmolbiophysbc.org/publications#asyncre_software_2015)

Installation & Usage
--------------------

It is recommended that the installation is performed in a personal python environment (`miniforge`, `miniconda`, conda`, or similar). AToM-OpenMM requires the `openmm`, `configobj` and `numpy` python modules. 


This version of AToM-OpenMM requires OpenMM 8.1.1 or later. This conda command installs the necessary requirements:
```
conda create -n atm8.1.2 -c conda-forge openmm=8.1.2 ambertools openmmforcefields configobj setproctitle r-base
conda activate atm8.1.2
```
`setproctitle` above is optional but useful to track the names of the processes started by AToM-OpenMM. The `ambertools` package is not an actual dependency but it is needed to set up some of the systems in the examples. `openmmforcefields` is used for force field parameter assignments using OpenFF. `r-base` with the `UWHAM R package` (see below) is required for free energy estimation. See [examples](examples/) for examples and tutorials.

Finally, install AToM-OpenMM:

- From the latest sources:
```
git clone git clone https://github.com/Gallicchio-Lab/AToM-OpenMM.git
pip install .
```

- From the latest release:
```
wget https://github.com/Gallicchio-Lab/AToM-OpenMM/archive/refs/tags/v8.1.2.tar.gz
tar zxvf v8.1.2.tar.gz
cd AToM-OpenMM-8.1.2
pip install .
```

And this will install the UWHAM R package:
```
Rscript -e 'install.packages("UWHAM", repos = "http://cran.us.r-project.org")' 
```

While we strive to develop and distribute high-quality and bug-free software, keep in mind that this is research software under heavy development. AToM-OpenMM is provided without any guarantees of correctness. Please report issues [here](https://github.com/Gallicchio-Lab/AToM-OpenMM/issues). We welcome contributions and pull requests.

Documentation
-------------

[AToM-OpenMM User Guide](https://www.compmolbiophysbc.org/atom-openmm)

See [examples](examples/) for examples and tutorials.

Licensing
---------

 This software is licensed under the terms of the [GNU General Public License](http://opensource.org/licenses/GPL-3.0). See [LICENSE](LICENSE)
