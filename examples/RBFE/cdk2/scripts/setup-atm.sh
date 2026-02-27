#!/bin/bash

#load settings
work_dir=$(pwd)
scripts_dir=${work_dir}/scripts
. ${scripts_dir}/setup-settings.sh || exit 1

#get the alignment atoms
#generates alignments.pkl file in the ligands folder
cd ${work_dir}/ligands || exit 1
python ${scripts_dir}/find_alignment_atoms.py --refligFile ${ref_ligand}.sdf --LIG1refatoms "${ref_ligand_alignment_atoms}" --alignmentsYAMLoutFile alignments.yaml  || exit 1

cd ${work_dir} || exit 1

#process each ligand pair
npairs=${#ligands[@]}
n=$(expr ${#ligands[@]} - 1)
for l in `seq 0 $n` ; do
    
    #return to main work directory for each ligand pair
    cd ${work_dir} || exit 1

    #ligands in ligand pair
    ligpair=${ligands[$l]}
    read -ra ligs <<< $ligpair
    lig1=${ligs[0]}
    lig2=${ligs[1]}

    echo "Processing ligand pair $lig1 $lig2 ..."
    
    #append to list of ligand pairs
    if [ -z "$ligpreppairs" ] ; then
	ligpreppairs="${lig1}-${lig2}"
    else
	ligpreppairs="${ligpreppairs} ${lig1}-${lig2}"
    fi

    #creates system in complexes folder
    jobname=${receptor}-${lig1}-${lig2}
    jobwdir=${work_dir}/complexes/${jobname} 
    mkdir -p ${jobwdir} || exit 1
    cd ${jobwdir}  || exit 1

    #slurm/bash batch script
    sed "s#<JOBNAME>#${jobname}#g ; s#<WORKDIR>#${jobwdir}# ; s#<SCRIPTS_DIR>#${scripts_dir}# ; s#<CONDADIR>#${CONDA_PREFIX}# ;  s#<CONDAENV>#${CONDA_DEFAULT_ENV}# ; s#<RCPT>#${receptor}# ; s#<LIG1>#${lig1}# ; s#<LIG2>#${lig2}# " < ${scripts_dir}/run_template.sh > ${jobwdir}/run.sh

    cp ${scripts_dir}/{vmd_template.in} ${jobwdir}/

done
