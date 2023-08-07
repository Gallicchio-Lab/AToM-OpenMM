AToM-OpenMM v8.0.0rc
====================

The Alchemical Transfer Method for OpenMM (AToM-OpenMM) is an extensible Python package for the estimation of absolute and relative binding free energies of molecular complexes. It implements the [Alchemical Transfer Method (ATM)](https://pubs.acs.org/doi/10.1021/acs.jcim.1c01129) with  asynchronous parallel replica exchange molecular dynamics with the [OpenMM](https://github.com/openmm) library. The AToM software can be deployed on workstations or cluster nodes with one or more GPUs.

This version of AToM uses the [ATMForce potential](https://github.com/openmm/openmm/pull/4110) in the latest [OpenMM sources](https://github.com/openmm/openmm). It requires the compilation of OpenMM from sources.

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

It is recommended that the installation is performed in a personal python environment (`conda`, `miniconda`, or similar). AToM-OpenMM requires the `openmm`, `configobj` and `numpy` python modules. 

This conda command installs the necessary requirements to run the examples:
```
conda create -n atm -c conda-forge ambertools openmmforcefields configobj setproctitle r-base
```
`setproctitle` above is optional but useful to track the names of the processes started by AToM-OpenMM. The `ambertools` package is not an actual dependency but it is needed to set up some of the systems in the examples. `openmmforcefields` is used for force field parameter assignments using OpenFF. `r-base` with the `UWHAM R package` (see below) is required for free energy estimation. See [examples](examples/) for examples and tutorials.

This version of AToM-OpenMM requires the latest development version of OpenMM which, as of this writing, has not yet been released on conda-forge. Follow [OpenMM's installation instructions from sources](http://docs.openmm.org/latest/userguide/library/02_compiling.html). Note that the `conda create` command above has likely installed a released version of OpenMM, which would need to be overwritten or bypassed. One approach is to install OpenMM from sources in a chosen location and bypass the version installed in conda. If OpenMM is installed under `$HOME/devel/openmm` and OpenMM's python libraries are installed under `$HOME/devel/openmm/lib/python3.7/site-packages`, the following settings should get you going:
```
export OPENMM_DIR=$HOME/devel/openmm
export OPENMM_PLUGIN_DIR=${OPENMM_DIR}/lib/plugins
export LD_LIBRARY_PATH=${OPENMM_DIR}/lib:${OPENMM_DIR}/lib/plugins:$LD_LIBRARY_PATH
export PYTHONPATH=${OPENMM_DIR}/lib/python3.7/site-packages:$PYTHONPATH
```

Finally, install AToM-OpenMM:
```
git clone https://github.com/Gallicchio-Lab/AToM-OpenMM.git
cd AToM-OpenMM
git checkout v8.0.0rc
python setup.py install
```

And this will install the UWHAM R package:
```
Rscript -e 'install.packages("UWHAM", repos = "http://cran.us.r-project.org")' 
```

While we strive to develop and distribute high-quality and bug-free software, keep in mind that this is research software under heavy development. AToM-OpenMM is provided without any guarantees of correctness. Please report issues [here](https://github.com/Gallicchio-Lab/AToM-OpenMM/issues). We welcome contributions and pull requests.

Documentation
-------------

[Under construction](https://www.compmolbiophysbc.org/atom-openmm)

Licensing
---------

 This software is licensed under the terms of the [GNU General Public License](http://opensource.org/licenses/GPL-3.0). See [LICENSE](LICENSE)
