# API Reference

This page documents the installed AToM-OpenMM command-line entry points. It is organized as one page so the sidebar stays uncluttered; use the page table of contents to jump to a command.

For workflow-level background, start with the [User Guide](../user-guide/index.md). For YAML configuration files, see [Configuration Files](../user-guide/configuration.md).

**Command summary**

| <span style="display: inline-block; min-width: 10rem;">Command</span> | Use it when you want to... | Main input | Main output |
| --- | --- | --- | --- |
| [`abfe_structprep`](#abfe_structprep) | Prepare an ABFE system for production. | ABFE YAML, JSON, or CNTL file. | Equilibrated `_0.xml` and `_0.pdb` files. |
| [`rbfe_structprep`](#rbfe_structprep) | Prepare an RBFE system for production. | RBFE YAML, JSON, or CNTL file. | Equilibrated `_0.xml` and `_0.pdb` files. |
| [`abfe_production`](#abfe_production) | Run ABFE asynchronous replica exchange. | ABFE YAML, JSON, or CNTL file. | Replica folders, trajectories, perturbation-energy data. |
| [`rbfe_production`](#rbfe_production) | Run RBFE asynchronous replica exchange. | RBFE YAML, JSON, or CNTL file. | Replica folders, trajectories, perturbation-energy data. |
| [`make_atm_system_from_pdb`](#make_atm_system_from_pdb) | Build an OpenMM system from a prepared system PDB and ligand SDF data. | System PDB plus ligand/cofactor SDF. | OpenMM System XML and topology PDB. |
| [`make_atm_system_from_rcpt_lig`](#make_atm_system_from_rcpt_lig) | Build and solvate a system from separate receptor and ligand files. | Receptor file, ligand file(s), displacement. | OpenMM System XML and topology PDB. |
| [`make_atm_system_from_amber`](#make_atm_system_from_amber) | Convert Amber topology/coordinates to OpenMM files. | Amber `prmtop` and `inpcrd`. | OpenMM System XML and topology PDB. |
| [`uwham`](#uwham) | Computes free energy estimates from production output. | Production run directory and job basename. | Free energy estimate, optional CSV files and plot. |

## `abfe_structprep`

Prepares an absolute binding free energy system for production. It reads the configuration, loads `<BASENAME>.pdb` and `<BASENAME>_sys.xml`, performs minimization, thermalization, NPT/NVT equilibration, lambda annealing, and final equilibration at the alchemical intermediate.

**Usage**

```bash
abfe_structprep system.yaml
```

The input file can be YAML, JSON, or CNTL. Newer workflows usually generate a per-system YAML file in the calculation directory.

**Expected inputs**

- A configuration file containing at least the job `BASENAME` and ATM settings.
- `<BASENAME>.pdb`, the topology/coordinate PDB for the system.
- `<BASENAME>_sys.xml`, the serialized OpenMM `System`.

**Key outputs**

| File pattern | Meaning |
| --- | --- |
| `<BASENAME>_min.xml`, `<BASENAME>_min.pdb` | Minimized system. |
| `<BASENAME>_therm.xml`, `<BASENAME>_therm.pdb` | Thermalized system. |
| `<BASENAME>_npt.xml`, `<BASENAME>_npt.pdb` | NPT-equilibrated system. |
| `<BASENAME>_equil.xml`, `<BASENAME>_equil.pdb` | NVT-equilibrated system before lambda annealing. |
| `<BASENAME>_mdlambda.xml`, `<BASENAME>_mdlambda.pdb`, `<BASENAME>_mdlambda.out` | Lambda-annealed state. |
| `<BASENAME>_0.xml`, `<BASENAME>_0.pdb`, `<BASENAME>_0.xtc` | Final prepared state used by production. |

**Notes**

- Run this before `abfe_production`.
- The command changes into the directory containing the input file, so relative paths should be written from that directory.
- By default, non-solvent atoms are temporarily restrained during minimization and thermalization. This behavior is controlled by `MINTHERM_RESTRAIN_SOLUTES`.

## `rbfe_structprep`

Prepares a relative binding free energy system for production. It follows the same broad preparation stages as `abfe_structprep`, but uses RBFE-specific system construction and ATM force setup.

**Usage**

```bash
rbfe_structprep system.yaml
```

**Expected inputs**

- A configuration file containing RBFE atom selections, alchemical schedule, displacement, restraints, and `BASENAME`.
- `<BASENAME>.pdb`, the topology/coordinate PDB for the paired system.
- `<BASENAME>_sys.xml`, the serialized OpenMM `System`.

In newer helper workflows, `run-atm.py` generates the per-system YAML file and fills in atom selections such as ligand atom lists, attachment atoms, receptor-frame atoms, alignment atoms, displacement, and exclusion-potential atom lists.

**Key outputs**

The output naming pattern matches `abfe_structprep`:

- `<BASENAME>_min.*`
- `<BASENAME>_therm.*`
- `<BASENAME>_npt.*`
- `<BASENAME>_equil.*`
- `<BASENAME>_mdlambda.*`
- `<BASENAME>_0.xml`, `<BASENAME>_0.pdb`, `<BASENAME>_0.xtc`

**Notes**

- Run this before `rbfe_production`.
- For small-molecule RBFE, alignment atom settings are commonly generated from ligand alignment inputs.
- For protein-peptide RBFE, atom selections can be generated from the mutation residue and peptide-backbone settings.

## `abfe_production`

Runs ABFE asynchronous replica exchange production. It reads the configuration, builds the ABFE replica-exchange job, dispatches OpenMM workers, writes trajectory and perturbation-energy data, and checkpoints the calculation.

**Usage**

```bash
abfe_production system.yaml
```

**Expected inputs**

- The same per-system configuration used for structure preparation.
- Prepared files from `abfe_structprep`, especially `<BASENAME>_0.xml`.
- A `nodefile` or suitable platform settings, depending on the workflow and hardware setup.

**Key outputs**

| Output | Meaning |
| --- | --- |
| `r0/`, `r1/`, ... | Replica directories. |
| `r*/<BASENAME>.out` | Perturbation-energy and state data used for analysis. |
| `r*/<BASENAME>.xtc` | Replica trajectory files. |
| `r*/<BASENAME>_ckpt.xml` | Replica checkpoints. |
| `<BASENAME>_stat.txt` | Replica-exchange status/statistics file. |

**Notes**

- The production command can restart from checkpoint files when they are present and valid.
- Sampling length is controlled by keys such as `WALL_TIME`, `MAX_SAMPLES`, `CYCLE_TIME`, and `CHECKPOINT_TIME`.
- Use `uwham` after production to recompute free energy estimates or write analysis CSV files.

## `rbfe_production`

Runs RBFE asynchronous replica exchange production. This is the RBFE counterpart to `abfe_production`; it uses RBFE system definitions and writes the same style of replica output.

**Usage**

```bash
rbfe_production system.yaml
```

**Expected inputs**

- The same per-system configuration used for RBFE structure preparation.
- Prepared files from `rbfe_structprep`, especially `<BASENAME>_0.xml`.
- A `nodefile` or suitable platform settings for the hardware being used.

**Key outputs**

| Output | Meaning |
| --- | --- |
| `r0/`, `r1/`, ... | Replica directories. |
| `r*/<BASENAME>.out` | Perturbation-energy and state data for both RBFE legs. |
| `r*/<BASENAME>.xtc` | Replica trajectory files. |
| `r*/<BASENAME>_ckpt.xml` | Replica checkpoints. |
| `<BASENAME>_stat.txt` | Replica-exchange status/statistics file. |

**Notes**

- For variable-displacement workflows, generated keywords such as `DISPLACEMENT`, `LIGOFFSET`, ligand atom lists, and exclusion atom lists are important. Inspect the generated YAML if behavior is unexpected.
- As with ABFE, production can continue from valid checkpoints.

## `make_atm_system_from_pdb`

Builds an OpenMM system from a prepared system PDB. Use this when the receptor, ligands, cofactors, solvent, and ions are already present in one PDB file, and nonstandard molecules are provided in an SDF file for force-field assignment.

**Usage**

```bash
make_atm_system_from_pdb \
  --systemPDBinFile system.pdb \
  --ligandsSDFFile ligands.sdf \
  --LIG1resid 123 \
  --systemXMLoutFile system_sys.xml \
  --systemPDBoutFile system.pdb
```

For RBFE-style inputs with two ligands, also provide:

```bash
  --LIG1refatoms "1 2 3" \
  --LIG2resid 124 \
  --LIG2refatoms "1 2 3"
```

**Important arguments**

| Argument | Required | Meaning |
| --- | --- | --- |
| `--systemPDBinFile` | Yes | Prepared PDB containing the full system. |
| `--ligandsSDFFile` | Yes | SDF containing all ligands and cofactors needing small-molecule parameters. |
| `--LIG1resid` | Yes | Residue ID of the first ligand in the PDB. |
| `--systemXMLoutFile` | Yes | Output OpenMM System XML file. |
| `--systemPDBoutFile` | Yes | Output topology/coordinate PDB file. |
| `--LIG1refatoms` | RBFE only | Reference atom IDs for ligand 1, starting from 1 within the ligand. |
| `--LIG2resid` | RBFE only | Residue ID of the second ligand. |
| `--LIG2refatoms` | RBFE only | Reference atom IDs for ligand 2. |
| `--proteinForceField` | No | Protein force field; default is `amber14-all.xml`. |
| `--solventForceField` | No | Solvent/ion force field; default is `amber14/tip3p.xml`. |
| `--ligandForceField` | No | Small-molecule force field; default is `openff-2.0.0`. |
| `--forcefieldJSONCachefile` | No | Cache file for generated ligand templates. |
| `--hmass` | No | Hydrogen mass; use values such as `1.5` for hydrogen mass repartitioning. |

**Outputs**

- The serialized OpenMM System XML requested by `--systemXMLoutFile`.
- The topology/coordinate PDB requested by `--systemPDBoutFile`.

## `make_atm_system_from_rcpt_lig`

Builds and solvates an ATM-ready OpenMM system from separate receptor and ligand input files. Use this for workflows that start from separate receptor and ligand structures rather than a fully assembled system PDB.

**Usage**

ABFE-style one-ligand setup:

```bash
make_atm_system_from_rcpt_lig \
  --receptorinFile receptor.pdb \
  --LIG1inFile ligand.sdf \
  --displacement "22.0 0.0 0.0" \
  --systemXMLoutFile system_sys.xml \
  --systemPDBoutFile system.pdb
```

RBFE-style two-ligand setup:

```bash
make_atm_system_from_rcpt_lig \
  --receptorinFile receptor.pdb \
  --LIG1inFile ligand1.sdf \
  --LIG2inFile ligand2.sdf \
  --displacement "22.0 0.0 0.0" \
  --systemXMLoutFile system_sys.xml \
  --systemPDBoutFile system.pdb
```

**Important arguments**

| Argument | Required | Meaning |
| --- | --- | --- |
| `--receptorinFile` | Yes | Receptor in PDB or SDF format. |
| `--displacement` | Yes | Displacement vector in Angstroms, written as a quoted string. |
| `--systemXMLoutFile` | Yes | Output OpenMM System XML file. |
| `--systemPDBoutFile` | Yes | Output topology/coordinate PDB file. |
| `--LIG1inFile` | Usually | First ligand in SDF or PDB format. |
| `--LIG2inFile` | RBFE only | Second ligand in SDF or PDB format. |
| `--cofactorsSDFFile` | No | SDF file containing receptor cofactors. |
| `--proteinForceField` | No | Protein force field. Can be supplied more than once. |
| `--solventForceField` | No | Solvent/ion force field. Can be supplied more than once. |
| `--ligandForceField` | No | Ligand force field: OpenFF, GAFF, or Espaloma-style names. |
| `--implicitSolvent` | No | Implicit solvent model such as `HCT`, `OBC2`, `GBN2`, or `vacuum`. |
| `--forcefieldJSONCachefile` | No | Cache file for ligand templates. |
| `--hmass` | No | Hydrogen mass. |
| `--ionicStrength` | No | Monatomic ion concentration for explicit solvent; default is `0.15`. |

**Outputs**

- A solvated or implicit-solvent OpenMM System XML.
- A matching PDB file containing the assembled system.

**Notes**

- If only `--LIG1inFile` is supplied, the command builds a one-ligand system.
- If `--LIG2inFile` is supplied, the second ligand is placed using the displacement vector.
- Deprecated aliases `--LIG1SDFinFile` and `--LIG2SDFinFile` still exist, but new workflows should use `--LIG1inFile` and `--LIG2inFile`.

## `make_atm_system_from_amber`

Converts an Amber topology and coordinate pair into the OpenMM files used by AToM-OpenMM. Use this when the system has already been prepared with AmberTools or another Amber-compatible workflow.

**Usage**

```bash
make_atm_system_from_amber \
  --AmberPrmtopinFile system.prmtop \
  --AmberInpcrdinFile system.inpcrd \
  --systemXMLoutFile system_sys.xml \
  --systemPDBoutFile system.pdb
```

**Arguments**

| Argument | Required | Meaning |
| --- | --- | --- |
| `--AmberPrmtopinFile` | Yes | Amber topology file. |
| `--AmberInpcrdinFile` | Yes | Amber coordinate file. |
| `--systemXMLoutFile` | Yes | Output OpenMM System XML file. |
| `--systemPDBoutFile` | Yes | Output topology/coordinate PDB file. |
| `--hmass` | No | Hydrogen mass; default is `1.0`. |
| `--nonbondedCutoff` | No | Nonbonded cutoff in nm; default is `0.9`. |
| `--switchDistance` | No | Switching distance in nm; default is `0.0`. |
| `--verbose` | No | Print additional output. |

**Outputs**

- The serialized OpenMM System XML requested by `--systemXMLoutFile`.
- The PDB requested by `--systemPDBoutFile`.

## `uwham`

Runs UWHAM analysis on completed AToM-OpenMM production data. It reads replica output files, estimates the free energy, and can write analysis CSV files and a quality-control plot.

**Usage**

```bash
uwham --jobname system
```

From outside the run directory:

```bash
uwham --jobname system --rundir /path/to/system
```

With explicit sample range and output files:

```bash
uwham \
  --jobname system \
  --mintimeid 20 \
  --leg1DataCSVoutFile leg1.csv \
  --leg2DataCSVoutFile leg2.csv \
  --plotOutFile qc.png
```

**Arguments**

| <span style="display: inline-block; min-width: 10rem;">Argument</span> | Required | Meaning |
| --- | --- | --- |
| `--jobname` | Yes | Basename of the production run. |
| `--rundir` | No | Directory containing the replica folders; defaults to the current directory. |
| `--mintimeid` | No | Only process samples after this time/sample ID. If omitted, the command discards the first third of samples. |
| `--maxtimeid` | No | Only process samples before this time/sample ID. |
| `--leg1DataCSVoutFile` | No | Write processed leg 1 samples and WHAM weights to CSV. |
| `--leg2DataCSVoutFile` | No | Write processed leg 2 samples and WHAM weights to CSV. |
| `--plotOutFile` | No | Write a PNG quality-assessment plot. |

**Expected inputs**

- Replica folders such as `r0/`, `r1/`, and so on.
- Production output files named `r*/<jobname>.out`.

**Output**

The command prints a free energy estimate like:

```text
system: DG = ... +/- ... DG(leg1) = ... DG(leg2) = ... kcal/mol, n_samples: ...
```

For ABFE-style workflows, helper scripts may label the final binding free energy as `DGb`. For RBFE-style workflows, helper scripts may label the relative binding free energy as `DDGb`.
