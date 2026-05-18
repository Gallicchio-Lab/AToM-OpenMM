Absolute Binding Free Energies of FKBP Ligands with Variable Displacement
-------------------------------------------------------------------------

This workflow calculates absolute binding free energies for FKBP ligand fragments from the work of [A Pan, H Xu, T Palpant and DE Shaw](http://dx.doi.org/10.1021/acs.jctc.7b00172). It mirrors the `examples/RBFE/cdk2` variable-displacement setup as closely as possible, but applies it to FKBP absolute binding free energies. Conceptually, it is RBFE variable-displacement reused for ABFE by replacing one ligand with a noninteracting ghost particle.

We assume that the latest release of OpenMM and the latest AToM-OpenMM packages and their dependencies are available in a `conda` environment (see the [README](../../../README.md), and that the `examples` folder is available under `$HOME/AToM-OpenMM/examples`. Adjust this pathname as needed.

The mapping is:

- `L1` is the physical ligand in the binding pocket.
- `L2` is a one-atom ghost residue with mass and no interactions.
- `run-atm.py` builds a standard one-ligand system, leaves the ligand in the binding pocket, appends the ghost at the displaced solvent position, and then feeds RBFE-style keywords into `rbfe_structprep` and `rbfe_production`.

The setup and launch flow intentionally follows the `cdk2` example:

```bash
cd $HOME/AToM-OpenMM/examples/ABFE/fkbp
bash ./scripts/setup-atm.sh
cd complexes

# run on a SLURM cluster
for i in fkbp-* ; do ( cd $i && sbatch ./run.sh ) ; done

# or run locally
# for i in fkbp-* ; do ( cd $i && bash ./run.sh > ${i}.log 2>&1 ) ; done
```

Collect the results once the simulations have completed:

```bash
cd $HOME/AToM-OpenMM/examples/ABFE/fkbp/complexes
for i in *; do tail -2 ${i}/${i}.log; done
```

The final log lines report the binding free energy estimate as `DGb` in kcal/mol, its estimated standard error, the ligand and ghost leg free energies, and the number of perturbation energy samples used after discarding the initial samples. These short tutorial calculations are configured for modest runtimes; significantly longer production runs, larger `MAX_SAMPLES` values, and correspondingly larger `WALL_TIME` and SLURM limits are needed for converged quantitative free energies.

Free energies and quality-control plots can also be regenerated from a completed calculation directory with the `uwham` application:

```bash
cd $HOME/AToM-OpenMM/examples/ABFE/fkbp/complexes/fkbp-dmso
uwham --jobname fkbp-dmso --plotOutFile fkbp-dmso.png
```

The quality-control plot compares the alchemical legs and the statistical weights used by UWHAM. Well-behaved calculations should show smooth free-energy trends and weights that are not dominated by a small number of samples.

After or while production is underway, view the alchemical trajectories with VMD using the generated structure and VMD script:

```bash
cd $HOME/AToM-OpenMM/examples/ABFE/fkbp/complexes/fkbp-dmso
vmd -f fkbp-dmso_0.pdb $(/bin/ls -v r*/*xtc) -e vmd.in
```

`scripts/setup-settings.sh` only lists the receptor basename and the ligand basenames. The initial displacement is determined automatically in `scripts/run-atm.py`, matching the variable-displacement behavior used by the RBFE workflow.

Notes:

- Positional restraints are disabled for this workflow.
- No alignment-restraint settings are needed in `scripts/defaults.yaml`.
- If you want a specific ligand attachment atom instead of the default heavy atom closest to the ligand centroid, pass `--ligandAttachAtomIndex` to `scripts/run-atm.py`.
