#!/bin/bash
#
#SBATCH -J <JOBNAME>
#SBATCH --partition=<PARTITION>
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --account=<ACCOUNTNO>
#SBATCH --no-requeue
#SBATCH -t 10:00:00

if [ -n "$CUDA_VISIBLE_DEVICES" ] ; then
    echo "localhost,0:$CUDA_VISIBLE_DEVICES,1,centos-OpenCL,,/tmp" > nodefile
fi
<SCRIPTS_DIR>/runopenmm <ASYNCRE_DIR>/rbfe_explicit.py <JOBNAME>_asyncre.cntl
