#!/bin/bash
#
#SBATCH -J atmbn-prep
#SBATCH --partition=gpu-shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --no-requeue
#SBATCH --account=<account>
#SBATCH -t 12:00:00

. ~/miniconda3/bin/activate atm
for pair in <LIGPAIRS> ; do
    jobname=<RECEPTOR>-$pair
    echo "Prepping $jobname"
    ( cd ${jobname} &&  python <ASYNCRE_DIR>/rbfe_structprep.py ${jobname}_asyncre.cntl )  || exit 1
done
