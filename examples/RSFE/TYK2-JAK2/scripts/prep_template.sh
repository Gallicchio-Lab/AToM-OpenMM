#!/bin/bash
#
#SBATCH -J atmbn-prep
#SBATCH --partition=any
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=2
#SBATCH --no-requeue
#SBATCH -t 12:00:00

. /nfs/sazimi-d/miniforge3/bin/activate atm8.2.0

for pair in <LIGPAIRS> ; do
    jobname=<RECEPTOR>-$pair
    echo "Prepping $jobname"
    ( cd ${jobname} &&  python <ASYNCRE_DIR>/rbfe_structprep.py ${jobname}_asyncre.cntl )  || exit 1
done
