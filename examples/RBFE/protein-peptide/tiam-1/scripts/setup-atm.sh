#!/bin/bash

#load settings
work_dir=$(pwd)
scripts_dir=${work_dir}/scripts
. ${scripts_dir}/setup-settings.sh || exit 1

workflow_mode=${workflow_mode:-ligand_rbfe}
ligand_structures_dir=${work_dir}/ligands
ligand_file_ext=sdf
alignment_file=${ligand_structures_dir}/alignments.yaml

case "${workflow_mode}" in
    ligand_rbfe)
        #generates alignments.yaml file in the ligands folder
        cd ${ligand_structures_dir} || exit 1
        python ${scripts_dir}/find_alignment_atoms.py --refligFile ${ref_ligand}.sdf --LIG1refatoms "${ref_ligand_alignment_atoms}" --alignmentsYAMLoutFile alignments.yaml || exit 1
        ;;
    protein_mutation)
        ligand_file_ext=pdb
        alignment_file=${ligand_structures_dir}/pair_alignments.yaml
        pair_args=()
        for ligpair in "${ligands[@]}" ; do
            read -ra ligs <<< "${ligpair}"
            if [ "${#ligs[@]}" -ne 3 ] ; then
                echo "Invalid protein_mutation pair: '${ligpair}'. Expected 'MUT1 MUT2 RESID'." >&2
                exit 1
            fi
            pair_args+=( --pair "${ligpair}" )
        done
        python ${scripts_dir}/find_pair_alignment_atoms.py \
            --jobPrefix "${receptor}" \
            --structuresDir "${ligand_structures_dir}" \
            --pairAlignmentsYAMLoutFile "${alignment_file}" \
            --fileExtension "${ligand_file_ext}" \
            "${pair_args[@]}" || exit 1
        ;;
    *)
        echo "Unsupported workflow_mode '${workflow_mode}'. Expected 'ligand_rbfe' or 'protein_mutation'." >&2
        exit 1
        ;;
esac

cd ${work_dir} || exit 1

#process each ligand pair
npairs=${#ligands[@]}
n=$(expr ${#ligands[@]} - 1)
for l in `seq 0 $n` ; do
    
    #return to main work directory for each ligand pair
    cd ${work_dir} || exit 1

    #ligands in ligand pair
    ligpair=${ligands[$l]}
    read -ra ligs <<< "${ligpair}"
    lig1=${ligs[0]}
    lig2=${ligs[1]}
    mutation_site=
    if [ "${workflow_mode}" = "protein_mutation" ] ; then
        mutation_site=${ligs[2]}
    fi

    if [ -n "${mutation_site}" ] ; then
        echo "Processing transformation ${lig1} ${lig2} at ${mutation_site} ..."
    else
        echo "Processing ligand pair ${lig1} ${lig2} ..."
    fi
    
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
    if [ "${workflow_mode}" = "protein_mutation" ] ; then
        alignment_cli="--pairAlignmentsYAMLinFile \${sdir}/../ligands/$(basename "${alignment_file}")"
    else
        alignment_cli="--alignmentsYAMLinFile \${sdir}/../ligands/$(basename "${alignment_file}")"
    fi
    sed "s#<JOBNAME>#${jobname}#g ; s#<WORKDIR>#${jobwdir}# ; s#<SCRIPTS_DIR>#${scripts_dir}# ; s#<CONDADIR>#${CONDA_PREFIX}# ; s#<CONDAENV>#${CONDA_DEFAULT_ENV}# ; s#<RCPT>#${receptor}# ; s#<LIG1>#${lig1}# ; s#<LIG2>#${lig2}# ; s#<WORKFLOW_MODE>#${workflow_mode}# ; s#<LIGFILEEXT>#${ligand_file_ext}# ; s#<ALIGNMENT_CLI>#${alignment_cli}# " < ${scripts_dir}/run_template.sh > ${jobwdir}/run.sh

    cp ${scripts_dir}/vmd_template.in ${jobwdir}/

done
