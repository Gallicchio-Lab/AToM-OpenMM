#!/bin/bash
#
#SBATCH -J <JOBNAME>
#SBATCH --partition=<PARTITION>
#SBATCH --qos=<QOS>
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --account=<ACCOUNTNO>
#SBATCH --no-requeue
#SBATCH -t 10:00:00


echo "localhost,0:0,1,centos-OpenCL,,/tmp" > nodefile

<SCRIPTS_DIR>/runopenmm <ASYNCRE_DIR>/abfe_explicit.py <JOBNAME>_asyncre.cntl
