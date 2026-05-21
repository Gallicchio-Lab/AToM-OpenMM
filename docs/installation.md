# Installation

This page covers the documentation environment only. It does not install the AToM-OpenMM runtime dependencies used for molecular simulations.

Create the docs environment from the repository root:

```bash
conda env create -f docs/environment.yml
```

Activate it:

```bash
conda activate atom-openmm-docs
```

The environment installs MkDocs and Material for MkDocs from `conda-forge`, so a global MkDocs installation is not required.

To check the site:

```bash
mkdocs build --strict
```
