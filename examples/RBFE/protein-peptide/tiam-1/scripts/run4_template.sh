#!/bin/bash
#
#SBATCH -J <JOBNAME>
#SBATCH --partition=any
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:4
#SBATCH --cpus-per-task=6
#SBATCH --no-requeue
#SBATCH -t 24:00:00

jobname=<JOBNAME>

. ${HOME}/miniforge3/bin/activate <CONDAENV>
echo "Running on $(hostname)"

if [ ! -f ${jobname}_0.xml ]; then
   rbfe_structprep ${jobname}_asyncre.cntl || exit 1
fi

echo "localhost,0:0,1,CUDA,,/tmp" > nodefile
echo "localhost,0:0,1,CUDA,,/tmp" >> nodefile
echo "localhost,0:0,1,CUDA,,/tmp" >> nodefile
echo "localhost,0:0,1,CUDA,,/tmp" >> nodefile
rbfe_production ${jobname}_asyncre.cntl
