#!/bin/bash
#
#SBATCH -J <JOBNAME>
#SBATCH --partition=any
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=3
#SBATCH --no-requeue
#SBATCH -t 26:00:00

jobname=<JOBNAME>

. ${HOME}/miniconda3/bin/activate atm8.4.0
echo "Running on $(hostname)"

if [ ! -f ${jobname}_0.xml ]; then
    abfe_structprep ${jobname}_asyncre.cntl || exit 1
fi

echo "localhost,0:0,1,CUDA,,/tmp" > nodefile
abfe_production ${jobname}_asyncre.cntl
