#!/bin/bash

#load settings
work_dir=$(pwd)
scripts_dir=${work_dir}/scripts
. ${scripts_dir}/setup-settings.sh || exit 1

cd ${work_dir} || exit 1

#process each ligand
n=$(expr ${#ligands[@]} - 1)
for l in `seq 0 $n` ; do

    cd ${work_dir} || exit 1

    lig=${ligands[$l]}

    echo "Processing ligand ${lig} ..."

    jobname=${receptor}-${lig}
    jobwdir=${work_dir}/complexes/${jobname}
    mkdir -p ${jobwdir} || exit 1
    cd ${jobwdir} || exit 1

    sed "s#<JOBNAME>#${jobname}#g ; s#<WORKDIR>#${jobwdir}# ; s#<SCRIPTS_DIR>#${scripts_dir}# ; s#<RCPT>#${receptor}# ; s#<LIG1>#${lig}# " < ${scripts_dir}/run_template.sh > ${jobwdir}/run.sh || exit 1
    chmod +x ${jobwdir}/run.sh || exit 1
done
