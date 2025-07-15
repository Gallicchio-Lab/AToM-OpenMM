#!/bin/bash

#load settings
work_dir=$(pwd)
scripts_dir=${work_dir}/scripts
source ${scripts_dir}/setup-settings.sh

cd ${work_dir} || exit 1

#process each ligand
nligands=${#ligands[@]}
nlig_m1=$(expr ${#ligands[@]} - 1)

cd ${work_dir} || exit 1
for l in `seq 0 ${nlig_m1}` ; do

    lig=${ligands[$l]}

    echo "Processing ligand  $lig..."

    #creates system in complexes folder
    jobname=${receptor}-${lig}
    mkdir -p ${work_dir}/complexes/${jobname} || exit 1
    
    cd ${work_dir}/complexes/${jobname}    || exit 1

    rcptpdb=${work_dir}/receptor/${receptor}.pdb
    ligsdf=${work_dir}/ligands/${lig}.sdf

    displs=""
    count=0
    for d in ${displacement[@]}; do
	displs="$displs $d"
    done
    make_atm_system_from_rcpt_lig --receptorinFile "${rcptpdb}" --LIG1SDFinFile "${ligsdf}" --displacement "${displs}" --systemXMLoutFile "${jobname}_sys.xml" --systemPDBoutFile "${jobname}.pdb" --forcefieldJSONCachefile "${work_dir}/ligands/ffdb.json" || exit 1

    #residue ligand name
    lig1resname=L1

    vres=${vsite_rcpt_residues[@]}

    echo "python $scripts_dir/create_cntlfile_from_template.py --systemPDBFile "${jobname}.pdb" --templatein ${work_dir}/scripts/asyncre_template.cntl --jobname "${jobname}" --displacement "${displs}" --lig1resname "${lig1resname}" --vsiteResidues "${vres}" "
    
    python $scripts_dir/create_cntlfile_from_template_abfe.py --systemPDBFile "${jobname}.pdb" --templatein ${work_dir}/scripts/asyncre_template.cntl --jobname "${jobname}" --displacement "${displs}" --lig1resname "${lig1resname}" --vsiteResidues "${vres}" --cntlfileout ${jobname}_asyncre.cntl || exit 1

    #copy slurm files, etc
    sed "s#<JOBNAME>#${jobname}#" < ${work_dir}/scripts/run_template.sh > ${work_dir}/complexes/${jobname}/run.sh
    cp ${work_dir}/scripts/analyze.sh ${work_dir}/scripts/uwham_analysis.R ${work_dir}/complexes/${jobname}/
    
done

#prepare prep script
ligs=${ligands[@]}
sed "s#<RECEPTOR>#${receptor}#;s#<LIGS>#${ligs}#g "< ${work_dir}/scripts/prep_template.sh > ${work_dir}/complexes/prep.sh

#prepare free energy calculation script
sed "s#<RECEPTOR>#${receptor}#;s#<LIGS>#${ligs}# " < ${work_dir}/scripts/free_energies_template.sh > ${work_dir}/complexes/free_energies.sh
