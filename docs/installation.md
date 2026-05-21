# Installation

Install AToM-OpenMM in a personal Python environment such as Miniforge, Miniconda, conda, or mamba.

## Runtime Environment

AToM-OpenMM requires OpenMM 8.4.0 or later, along with Python dependencies such as `configobj` and `numpy`. A typical conda-forge environment can be created with:

```bash
mamba create -n atm8.4.0 -c conda-forge "openmm>=8.4" ambertools openmmforcefields configobj setproctitle r-base espaloma
mamba activate atm8.4.0
```

`setproctitle` is optional, but useful for tracking AToM-OpenMM processes. The `ambertools`, `openmmforcefields`, and `espaloma` packages are used to prepare molecular systems and assign force field parameters. `r-base` is used for installing the UWHAM R package used in free energy analysis.

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
