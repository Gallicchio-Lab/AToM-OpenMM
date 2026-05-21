# Installation

It is recommended that the installation is performed in a personal python environment (`miniforge`, `miniconda`, `conda`, or similar). AToM-OpenMM requires the `openmm`, `configobj` and `numpy` python modules. 

## Runtime Environment

This version of AToM-OpenMM requires OpenMM 8.4.0 or later. This conda command installs the necessary requirements:
```bash
mamba create -n atm8.4.0 -c conda-forge 'openmm>=8.4' ambertools openmmforcefields configobj setproctitle r-base espaloma
mamba activate atm8.4.0
```
`setproctitle` above is optional but useful to track the names of the processes started by AToM-OpenMM. The `ambertools`, `openmmforcefields`, and `espaloma` packages are not actual dependencies; they are used to setup the molecular systems. `openmmforcefields`, in particular, is used to assign force field parameters using OpenFF, GAFF, or `espaloma`. [`espaloma`](https://github.com/choderalab/espaloma) is a machine-learning system by the Chodera lab to assign force field parameters.  `r-base` with the `UWHAM R package` (see below) is required for free energy estimation. See [examples](examples/) for examples and tutorials.
## Install AToM-OpenMM

Install the latest release:

```bash
pip install atom-openmm
```

Or install from source:

```bash
git clone https://github.com/Gallicchio-Lab/AToM-OpenMM.git
cd AToM-OpenMM
pip install .
```

Install UWHAM for free energy estimation:

```bash
Rscript -e 'install.packages("UWHAM", repos = "http://cran.us.r-project.org")'
```
