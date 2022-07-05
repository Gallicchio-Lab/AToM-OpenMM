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

for lig in <LIGS> ; do
    jobname=<RECEPTOR>-${lig}
    echo "Prepping $jobname"
    ( cd ${jobname} &&  python ${jobname}_mintherm.py && python ${jobname}_equil.py && python ${jobname}_mdlambda.py  )  || exit 1
done
