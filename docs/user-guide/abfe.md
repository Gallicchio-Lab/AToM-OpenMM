# Absolute Binding Free Energy

Absolute binding free energy (ABFE) workflows estimate the binding free energy of one ligand or molecular partner against a target receptor. The current YAML-based ABFE workflow in this repository is `examples/ABFE/fkbp`, which adapts the newer variable-displacement RBFE machinery to an absolute binding calculation.

In this workflow, `L1` is the physical ligand in the binding site and `L2` is a one-atom noninteracting ghost particle placed in the displaced solvent position. The setup script builds a standard receptor-ligand system, patches in the ghost particle, writes the final per-system YAML file, runs structure preparation, runs production sampling, and analyzes the two alchemical legs with UWHAM.

## When to Use ABFE

Use ABFE when each ligand can be treated independently against the same receptor or host. This is useful when ligands are structurally diverse, or when you want an absolute binding estimate for a single complex.

Use RBFE instead when the main scientific question is a difference between two related ligands, peptides, or variants. RBFE is usually the more direct choice for ranking a congeneric series.

## Input Layout

The newer FKBP ABFE example uses this layout:

```text
examples/ABFE/fkbp/
  receptor/
    fkbp.pdb
  ligands/
    ligand-name.sdf
  scripts/
    setup-settings.sh
    defaults.yaml
    setup-atm.sh
    run-atm.py
    run_template.sh
  complexes/
    fkbp-ligand-name/
      run.sh
```

Before adapting the workflow, prepare the receptor and ligand structures carefully. The receptor PDB should already be chemically complete for the intended force field, and ligand structures should have sensible protonation, chemistry, and binding poses.

## Setup Files

`scripts/setup-settings.sh` defines the receptor basename and the ligand basenames to process. In the FKBP example it is intentionally short: one receptor name and a list of ligands.

`scripts/defaults.yaml` holds the shared simulation settings: force-field choices, the alchemical schedule, restraint behavior, and production length. See [Configuration Files](configuration.md) for the keyword-level reference.

The ABFE workflow does not need ligand alignment-restraint settings because the second alchemical partner is a one-atom ghost.

## Running the Workflow

From the example directory, generate calculation folders and launch each generated `run.sh`:

```bash
cd $HOME/AToM-OpenMM/examples/ABFE/fkbp
bash ./scripts/setup-atm.sh
cd complexes

# On a SLURM cluster
for i in fkbp-*; do ( cd "$i" && sbatch ./run.sh ); done

# Or without SLURM
# for i in fkbp-*; do ( cd "$i" && bash ./run.sh > "${i}.log" 2>&1 ); done
```

Each `run.sh` calls `scripts/run-atm.py` with the shared defaults, receptor file, ligand file, job name, and output plot name. The script chooses a displacement vector if one is not supplied, adds the ghost particle, computes atom selections, writes `<jobname>.yaml`, runs `rbfe_structprep`, runs `rbfe_production`, and then performs UWHAM analysis if production samples are present.

## Outputs to Check

Inside each `complexes/<jobname>/` directory, the important files are:

| Output | Use |
| --- | --- |
| `<jobname>.yaml` | Final merged options used for this prepared calculation. |
| `<jobname>.pdb`, `<jobname>_sys.xml` | Prepared system with the physical ligand and ghost particle. |
| `<jobname>_0.xml`, `<jobname>_0.pdb` | Prepared state used to start production. |
| `r*/<jobname>.out` | Perturbation-energy samples from each replica. |
| `<jobname>.log` | Runtime log and final free-energy summary. |
| `<jobname>.png` | UWHAM quality-control plot when requested. |
| `vmd.in` | VMD helper file generated from the template when available. |

The final log lines report `DGb` in kcal/mol, an estimated uncertainty, the ligand and ghost leg free energies, and the number of perturbation-energy samples used after discarding initial samples.

Analysis can be regenerated later from a completed calculation directory:

```bash
cd $HOME/AToM-OpenMM/examples/ABFE/fkbp/complexes/fkbp-dmso
uwham --jobname fkbp-dmso --plotOutFile fkbp-dmso.png
```

## Adapting to a New ABFE System

Start by copying the FKBP workflow pattern rather than editing generated calculation folders by hand. Update the receptor and ligand inputs, edit `setup-settings.sh`, and then tune `defaults.yaml` for the target system.

For a new ABFE calculation, pay special attention to the receptor chain definitions, ligand force-field choice, binding-site restraint behavior, alchemical schedule, soft-core settings, production length, and scheduler time limit. The tutorial defaults are intentionally short and should be increased for quantitative work.

If the default ligand attachment atom is not appropriate, pass `--ligandAttachAtomIndex` through the run command to choose a specific 1-based ligand atom index.

Quantitative ABFE estimates require longer sampling, convergence checks across the alchemical legs, and careful inspection of the quality-control plot and replica outputs.

## Older Fixed-Displacement ABFE

Some tutorials, such as `examples/ABFE/temoa-g1`, use the older `.cntl` workflow instead of the newer YAML workflow. This style is still supported and is useful when reproducing older examples or when starting from an already prepared Amber topology and coordinate file.

The older ABFE workflow looks like:

```bash
cd $HOME/AToM-OpenMM/examples/ABFE/temoa-g1
make_atm_system_from_amber \
    --AmberPrmtopinFile temoa-g1.prmtop \
    --AmberInpcrdinFile temoa-g1.inpcrd \
    --systemXMLoutFile temoa-g1_sys.xml \
    --systemPDBoutFile temoa-g1.pdb
abfe_structprep temoa-g1_asyncre.cntl
abfe_production temoa-g1_asyncre.cntl
./analyze.sh 20
```

In this style, the control file contains the final runtime settings directly. Instead of generating atom selections into a YAML file, users provide the atom selections, displacement information, restraints, and alchemical schedule in the control file. Execution devices are listed in a `nodefile`.

The fixed-displacement style places the ligand at a specified bulk-solvent displacement vector rather than using the newer variable-displacement machinery. For new ABFE work, the YAML ghost-ligand workflow is usually easier to adapt, but the `.cntl` workflow remains a compact reference for manually prepared systems.
