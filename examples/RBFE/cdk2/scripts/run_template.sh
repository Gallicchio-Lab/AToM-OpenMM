#!/bin/bash
#
#SBATCH -J <JOBNAME>
##SBATCH --partition=debug
#SBATCH --partition=project
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
##SBATCH --qos=restrained
#SBATCH --qos=maxjobs
#SBATCH --no-requeue
#SBATCH --account=insite
#SBATCH -t 24:00:00

. ~/miniconda3/bin/activate atm
echo "Running on $(hostname)"

if [ ! -f <JOBNAME>_0.xml ]; then
   python <JOBNAME>_mintherm.py && python <JOBNAME>_mdlambda.py || exit 1
fi

echo "localhost,0:0,1,CUDA,,/tmp" > nodefile
python <ASYNCRE_DIR>/rbfe_explicit.py <JOBNAME>_asyncre.cntl
