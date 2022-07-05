#!/bin/bash
#
#SBATCH -J <RECEPTOR>
#SBATCH --partition=gpu-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus=1
#SBATCH --account=<ACCOUNTNO>
#SBATCH --no-requeue
#SBATCH -t 02:15:00

for pair in <LIGPAIRS> ; do
    jobname=<RECEPTOR>-$pair
    echo "Prepping $jobname"
    ( cd ${jobname} &&  python ${jobname}_mintherm.py && python ${jobname}_equil.py  )  || exit 1
done
