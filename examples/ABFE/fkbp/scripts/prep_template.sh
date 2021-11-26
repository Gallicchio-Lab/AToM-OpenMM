#!/bin/bash
#SBATCH -J <JOBNAME>
#SBATCH -o logs/job-%A.%a.log
#SBATCH -N 1 # 1 node
#SBATCH -n 1 # 1 task
#SBATCH -c 2
#SBATCH --gres=gpu:1
#SBATCH --qos=maxjobs
#SBATCH --time=3:00:00 #runtime limit of 10 hours
#SBATCH -p project
#SBATCH -A insite
#SBATCH --reservation=advsim
#SBATCH --array=0-<NLIGANDS>

scripts_dir=<SCRIPTS_DIR>
complexdirs=($(find ./* -maxdepth 0 -type d | sort -V))
numjobs=${#complexdirs[@]}
complex_name=${complexdirs[SLURM_ARRAY_TASK_ID]}

echo "Prepping ${complex_name}"
cd ${complex_name} || exit 1

${scripts_dir}/runopenmm ${complex_name}_mintherm.py || exit 1
${scripts_dir}/runopenmm ${complex_name}_equil.py || exit 1
${scripts_dir}/runopenmm ${complex_name}_mdlambda.py || exit

cp ${complex_name}_0.xml ./leg1/${complex_name}_0.xml || exit 1
cp ${complex_name}_0_displaced.xml ./leg2/${complex_name}_0.xml || exit 1

for leg in leg1 leg2 ; do
    cp ${complex_name}.prmtop ./${leg}/ || exit 1
    cp ${complex_name}}.inpcrd ./${leg} || exit 1

done
