#!/bin/bash
#
#SBATCH -J <JOBNAME>
#SBATCH --partition=gpu-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --account=<ACCOUNTNO>
#SBATCH --no-requeue
#SBATCH -t 02:15:00

echo "localhost,0:0,1,OpenCL,,/tmp" > nodefile
../../scripts/runopenmm rbfe_explicit.py <JOBNAME>_asyncre.cntl
