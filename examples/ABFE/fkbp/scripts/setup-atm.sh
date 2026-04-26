#!/bin/bash

#load settings
work_dir=$(pwd)
scripts_dir=${work_dir}/scripts
. ${scripts_dir}/setup-settings.sh || exit 1

cd ${work_dir} || exit 1

#process each ligand
n=$(expr ${#ligands[@]} - 1)
for l in `seq 0 $n` ; do

    #return to main work directory for each ligand pair
    cd ${work_dir} || exit 1

    lig=${ligands[$l]}

    echo "Processing ligand ${lig} ..."

    #create system in complexes folder
    jobname=${receptor}-${lig}
    jobwdir=${work_dir}/complexes/${jobname}
    mkdir -p ${jobwdir} || exit 1
    cd ${jobwdir} || exit 1

    #slurm/bash batch script
    sed "s#<JOBNAME>#${jobname}#g ; s#<WORKDIR>#${jobwdir}# ; s#<SCRIPTS_DIR>#${scripts_dir}# ; s#<CONDADIR>#${CONDA_PREFIX}# ; s#<CONDAENV>#${CONDA_DEFAULT_ENV}# ; s#<RCPT>#${receptor}# ; s#<LIG1>#${lig}# " < ${scripts_dir}/run_template.sh > ${jobwdir}/run.sh || exit 1
    
    cp ${scripts_dir}/vmd_template.in ${jobwdir}/
done
