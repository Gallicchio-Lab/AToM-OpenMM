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

. $HOME/miniforge3/bin/activate
echo "localhost,0:0,1,CUDA,,/tmp" > nodefile
abfe_production <JOBNAME>_asyncre.cntl
