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
    ( cd ${jobname} &&  ../../scripts/runopenmm ${jobname}_mintherm.py && ../../scripts/runopenmm ${jobname}_equil.py && ../../scripts/runopenmm ${jobname}_mdlambda.py  )  || exit 1
done
