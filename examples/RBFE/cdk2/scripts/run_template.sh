#!/bin/bash
#
#SBATCH -J <JOBNAME>
#SBATCH --partition=gpu-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --account=TG-MCB150001
#SBATCH --no-requeue
#SBATCH -t 10:00:00

. ~/miniconda3/bin/activate atm
echo "Running on $(hostname)"

if [ ! -f <JOBNAME>_0.xml ]; then
   python <JOBNAME>_mintherm.py && python <JOBNAME>_mdlambda.py || exit 1
fi

echo "localhost,0:0,1,CUDA,,/tmp" > nodefile
python <ASYNCRE_DIR>/rbfe_explicit.py <JOBNAME>_asyncre.cntl
