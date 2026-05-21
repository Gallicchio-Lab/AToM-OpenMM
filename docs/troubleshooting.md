# Troubleshooting

## Environment Problems

AToM-OpenMM should be installed in an isolated conda or mamba environment. If imports fail, first confirm that the environment is active and that OpenMM is at least version 8.4.0:

```bash
python -c "import openmm; print(openmm.version.version)"
```

## OpenMM or GPU Issues

If OpenMM cannot find a GPU platform, check the CUDA drivers, the OpenMM installation, and the platform selected by your input files. On shared clusters, make sure the job has actually been assigned a GPU before starting production.

## UWHAM Is Missing

Free energy analysis uses the UWHAM R package. If analysis fails because UWHAM is unavailable, install it in the active environment:

```bash
Rscript -e 'install.packages("UWHAM", repos = "http://cran.us.r-project.org")'
```

## Example Paths

Many examples assume that the repository is available under `$HOME/AToM-OpenMM`. If your checkout lives elsewhere, adjust the paths in the tutorial commands and run scripts before submitting jobs.

## Reporting Issues

When reporting a problem, include the command you ran, the input files or configuration section involved, the relevant log output, and the OpenMM and AToM-OpenMM versions.
