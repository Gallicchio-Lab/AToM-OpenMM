#!/bin/bash
#
#SBATCH -J atmbn-prep
##SBATCH --partition=debug
#SBATCH --partition=project
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
##SBATCH --qos=restrained
#SBATCH --qos=maxjobs
#SBATCH --no-requeue
#SBATCH --account=insite
#SBATCH -t 12:00:00
. ~/miniconda3/bin/activate atm

for pair in <LIGPAIRS> ; do
    jobname=<RECEPTOR>-$pair
    echo "Prepping $jobname"
    ( cd ${jobname} &&  python ${jobname}_mintherm.py && python ${jobname}_mdlambda.py  )  || exit 1
done
