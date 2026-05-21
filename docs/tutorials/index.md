# Tutorials

Tutorials are the best place to start after reading the installation notes.

If you do not want to install AToM-OpenMM locally, start with the Colab notebooks. They provide guided examples in a browser-based environment and are useful for learning the workflow before setting up a workstation or cluster environment.

## Colab Notebooks

- [ABFE notebook](https://colab.research.google.com/github/Gallicchio-Lab/AToM-OpenMM/blob/master/example-notebooks/atom_openmm_abfe.ipynb)
- [RBFE notebook](https://colab.research.google.com/github/Gallicchio-Lab/AToM-OpenMM/blob/master/example-notebooks/atom_openmm_rbfe.ipynb)
- [Protein-peptide RBFE notebook](https://colab.research.google.com/github/Gallicchio-Lab/AToM-OpenMM/blob/master/example-notebooks/atom_openmm_protein_peptide_rbfe.ipynb)

## Local Tutorials

The repository tutorials are more detailed and are intended for local or cluster workflows. They are grouped by calculation type below.

### ABFE

ABFE tutorials estimate the binding free energy of one ligand or molecular partner against a receptor or host. The most up-to-date local ABFE workflow is the FKBP fragments example, which uses the newer YAML setup and ghost-ligand pattern described in the user guide.

- [ABFE FKBP fragments](https://github.com/Gallicchio-Lab/AToM-OpenMM/tree/master/examples/ABFE/fkbp): absolute binding free energies for FKBP ligand fragment complexes using the newer YAML workflow.
- [ABFE TEMOA-G1](https://github.com/Gallicchio-Lab/AToM-OpenMM/tree/master/examples/ABFE/temoa-g1): absolute binding free energy between the TEMOA host and the G1 guest from the SAMPL8 GDCC challenge.

### RBFE

RBFE tutorials estimate free energy differences between related ligands, guests, peptides, or molecular variants. For current local workflows, start with CDK2 for small-molecule ligand pairs and the TIAM1 protein-peptide example for mutation-like peptide variants.

- [RBFE CDK2 inhibitors](https://github.com/Gallicchio-Lab/AToM-OpenMM/tree/master/examples/RBFE/cdk2): relative binding free energies for CDK2 inhibitors using the newer YAML workflow and generated ligand alignment atoms.
- [Protein-peptide RBFE TIAM1 variants](https://github.com/Gallicchio-Lab/AToM-OpenMM/tree/master/examples/RBFE/protein-peptide/tiam1): relative binding free energies for Syndecan-1-derived peptide variants binding to the TIAM1 PDZ domain.
- [RBFE TEMOA G1-G4](https://github.com/Gallicchio-Lab/AToM-OpenMM/tree/master/examples/RBFE/temoa-g1-g4): relative binding free energy between the G1 and G4 guests to the TEMOA host.
- [RBFE ER alpha ligands](https://github.com/Gallicchio-Lab/AToM-OpenMM/tree/master/examples/RBFE/eralpha): relative binding free energies for ligands binding to the ER alpha protein receptor.

Each local tutorial assumes AToM-OpenMM and its runtime dependencies are installed. See [Installation](../installation.md) before running the repository examples.
