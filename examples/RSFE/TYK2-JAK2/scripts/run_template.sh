#!/bin/bash
#
#SBATCH -J <JOBNAME>
#SBATCH --partition=any
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --no-requeue
#SBATCH -t 24:00:00

jobname=<JOBNAME>

. /nfs/sazimi-d/miniforge3/bin/activate atm8.2.0

echo "Running on $(hostname)"

if [ ! -f ${jobname}_0.xml ]; then
   python <ASYNCRE_DIR>/rbfe_structprep.py ${jobname}_asyncre.cntl || exit 1
fi

echo "localhost,0:0,1,CUDA,,/tmp" > nodefile
python <ASYNCRE_DIR>/rbfe_explicit.py ${jobname}_asyncre.cntl
