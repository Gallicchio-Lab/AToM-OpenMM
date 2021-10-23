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

if [ -n "$CUDA_VISIBLE_DEVICES" ] ; then
    echo "localhost,0:$CUDA_VISIBLE_DEVICES,1,centos-OpenCL,,/tmp" > nodefile
fi
./runopenmm rbfe_explicit.py <JOBNAME>_asyncre.cntl
