Absolute Binding Free Energies of FKBP Ligands with Variable Displacement
-------------------------------------------------------------------------

This workflow mirrors the maintained [`examples/RBFE/cdk2`](../../RBFE/cdk2) variable-displacement setup as closely as possible, but applies it to FKBP absolute binding free energies. Conceptually, it is RBFE variable-displacement reused for ABFE by replacing one ligand with a noninteracting ghost particle.

The mapping is:

- `L1` is a one-atom ghost residue with mass and no interactions.
- `L2` is the physical ligand.
- `run-atm.py` builds a standard one-ligand system, translates the ligand to the displaced solvent site, appends the ghost at the original binding-site anchor position, and then feeds RBFE-style keywords into `rbfe_structprep` and `rbfe_production`.

The setup and launch flow intentionally follows the `cdk2` example:

```bash
cd $HOME/AToM-OpenMM/examples/ABFE/fkbp-displ
bash ./scripts/setup-atm.sh
cd complexes

# run on a SLURM cluster
for i in fkbp-* ; do ( cd $i && sbatch ./run.sh ) ; done

# or run locally
# for i in fkbp-* ; do ( cd $i && bash ./run.sh > ${i}.log 2>&1 ) ; done
```

`scripts/setup-settings.sh` only lists the receptor basename and the ligand basenames. The initial displacement is determined automatically in `scripts/run-atm.py`, matching the variable-displacement behavior used by the RBFE workflow.

Notes:

- Positional restraints are disabled for this workflow.
- No alignment-restraint settings are needed in `scripts/defaults.yaml`.
- The receptor and ligand inputs are copied from [`examples/ABFE/fkbp`](../fkbp).
- If you want a specific ligand attachment atom instead of the default first heavy atom, pass `--ligandAttachAtomIndex` to `scripts/run-atm.py`.
