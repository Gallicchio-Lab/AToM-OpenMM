#!/bin/bash
#
#SBATCH -J <JOBNAME>
#SBATCH --partition=<PARTITION>
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --account=<ACCOUNTNO>
#SBATCH --no-requeue
#SBATCH -t 02:15:00

python <ASYNCRE_DIR>/rbfe_explicit.py <JOBNAME>_asyncre.cntl
