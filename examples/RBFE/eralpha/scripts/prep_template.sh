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

for pair in <LEG1PAIRS> ; do
    mols_t=$(echo $pair | tr '-' ' ' )
    read -ra mols <<< $mols_t
    leg1=${mols[0]}-${mols[1]}
    leg2=${mols[1]}-${mols[0]}
    jobname1=<RECEPTOR>-$leg1
    jobname2=<RECEPTOR>-$leg2
    echo "Prepping $jobname1"
    cd ${jobname1} || exit 1
    ../../scripts/runopenmm ${jobname1}_mintherm.py || exit 1
    ../../scripts/runopenmm ${jobname1}_equil.py || exit 1
    cp ${jobname1}_0_displaced.xml ../${jobname2}/${jobname2}_0.xml || exit 1
    cp ${jobname1}.prmtop ../${jobname2}/${jobname2}.prmtop || exit 1
    cp ${jobname1}.inpcrd ../${jobname2}/${jobname2}.inpcrd || exit 1
    cd ..
done
