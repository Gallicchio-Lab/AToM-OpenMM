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

. $HOME/miniforge3/bin/activate
for lig in <LIGS> ; do
    jobname=<RECEPTOR>-${lig}
    echo "Prepping $jobname"
    ( cd ${jobname} &&  abfe_structprep ${jobname}_asyncre.cntl  )  || exit 1
done
