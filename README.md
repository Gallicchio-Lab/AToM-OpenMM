AToM-OpenMM
==============

The Alchemical Transfer Method for OpenMM (AToM-OpenMM) is an extensible Python package for the estimation of absolute and relative binding free energies of molecular complexes. It implements the [Alchemical Transfer Method (ATM)](https://pubs.acs.org/doi/10.1021/acs.jcim.1c01129) with  asynchronous parallel replica exchange molecular dynamics with the [OpenMM](https://github.com/openmm) library. The AToM software can be deployed on workstations or cluster nodes with one or more GPUs.

AToM uses the [ATMetaForce plugin](https://github.com/Gallicchio-Lab/openmm-atmmetaforce-plugin) for OpenMM.

Credits
-------

This software is developed and maintained by the [Emilio Gallicchio's lab](http://www.compmolbiophysbc.org) with support from current and past grants from the National Science Foundation (ACI 1440665 and CHE 1750511).

Authors:

Emilio Gallicchio <egallicchio@brooklyn.cuny.edu>

Baofeng Zhang BZhang@brooklyn.cuny.edu

Rajat Pal <rajatfor2014@gmail.com>

The asynchronous replica exchange method was first implemented in the [AsyncRE](https://github.com/ComputationalBiophysicsCollaborative/AsyncRE) package for the IMPACT program.

Citations
---------

Please [cite us](http://www.compmolbiophysbc.org/publications) if you use this software in your research:

- [Relative Binding Free Energy Calculations for Ligands with Diverse Scaffolds with the Alchemical Transfer Method](https://pubs.acs.org/doi/10.1021/acs.jcim.1c01129)

- [Alchemical Transfer Approach to Absolute Binding Free Energy Estimation](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00266)

- [Asynchronous Replica Exchange Software for Grid and Heterogeneous Computing](http://www.compmolbiophysbc.org/publications#asyncre_software_2015)

Installation & Usage
--------------------

It is recommended that the installation is performed in a personal python environment (`conda`, `miniconda`, or similar). AToM requires the `openmm`, `configobj` and `numpy` python modules. 

```
conda create -n atm -c conda-forge ambertools openmm openmm-atmmetaforce-plugin configobj setproctitle r-base
git clone https://github.com/Gallicchio-Lab/AToM-OpenMM.git
cd AToM-OpenMM
python setup.py install
Rscript -e 'install.packages("UWHAM", repos = "http://cran.us.r-project.org")' 
```

For NNP, install `openmm-ml` and `nnpops`:
```
conda install -c conda-forge openmm-ml nnpops
```

`setproctitle` above is optional but useful to track the names of the processes started by AToM-OpenMM. The `ambertools` package is not an actual dependency but it is needed to set up some of the systems in the examples. `r-base`, and the `UWHAM R package` is required for free energy estimation. See [examples](examples/) for examples and tutorials.

While we strive to develop and distribute high-quality and bug-free software, keep in mind that this is research software under heavy development. AToM is provided without any guarantees of correctness. Please report issues [here](https://github.com/Gallicchio-Lab/AToM-OpenMM/issues). We welcome contributions and pull requests.

Documentation
-------------

[Under construction](https://www.compmolbiophysbc.org/atom-openmm)

Licensing
---------

 This software is licensed under the terms of the [GNU General Public License](http://opensource.org/licenses/GPL-3.0). See [LICENSE](LICENSE)
