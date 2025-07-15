#!/bin/bash

#load settings
work_dir=$(pwd)
scripts_dir=${work_dir}/scripts
. ${scripts_dir}/setup-settings.sh

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
    mkdir -p ${work_dir}/complexes/${jobname} || exit 1
    cd ${work_dir}/complexes/${jobname}    || exit 1

    #files of receptor and ligands
    rcptpdb=${work_dir}/receptor/${receptor}.pdb
    lig1sdf=${work_dir}/ligands/${lig1}.sdf
    lig2sdf=${work_dir}/ligands/${lig2}.sdf
    
    #creates simulation system
    displs=""
    for d in ${displacement[@]}; do
	displs="$displs $d"
    done
    echo "Displacement vector: $displs"
    make_atm_system_from_rcpt_lig --receptorinFile "${rcptpdb}" --LIG1SDFinFile "${lig1sdf}" --LIG2SDFinFile "${lig2sdf}" --displacement "${displs}" --systemXMLoutFile "${jobname}_sys.xml" --systemPDBoutFile "${jobname}.pdb" --hmass 1.5 --forcefieldJSONCachefile "${work_dir}/ligands/ffdb.json" || exit 1
    #residue ligand names
    lig1resname=L1
    lig2resname=L2

    #ligand alignment reference atoms
    refpair=${ref_atoms[$l]}
    read -ra ref <<< $refpair
    ref_atoms1a="${ref[0]}"
    ref_atoms2a="${ref[1]}"

    #vsite residues
    vres=${vsite_rcpt_residues[@]}
    
    echo "python $scripts_dir/create_cntlfile_from_template_rbfe.py --systemPDBFile "${jobname}.pdb" --templatein ${work_dir}/scripts/asyncre_template.cntl --jobname "${jobname}" --displacement "${displs}" --lig1resname "${lig1resname}" --lig1refatoms "${ref_atoms1a}" --lig2resname "${lig2resname}" --lig2refatoms "${ref_atoms2a}"  --vsiteResidues "${vres}" --cntlfileout ${jobname}_asyncre.cntl "
    
    python $scripts_dir/create_cntlfile_from_template_rbfe.py --systemPDBFile "${jobname}.pdb" --templatein ${work_dir}/scripts/asyncre_template.cntl --jobname "${jobname}" --displacement "${displs}" --lig1resname "${lig1resname}" --lig1refatoms "${ref_atoms1a}" --lig2resname "${lig2resname}" --lig2refatoms "${ref_atoms2a}"  --vsiteResidues "${vres}" --cntlfileout ${jobname}_asyncre.cntl || exit 1

    sed "s#<JOBNAME>#${jobname}#g" < ${work_dir}/scripts/run_template.sh > ${work_dir}/complexes/${jobname}/run.sh

    cp ${work_dir}/scripts/analyze.sh ${work_dir}/scripts/uwham_analysis.R ${work_dir}/complexes/${jobname}/
    
done

#prepare prep script
sed "s#<RECEPTOR>#${receptor}# ; s#<LIGPAIRS>#${ligpreppairs}#" < ${work_dir}/scripts/prep_template.sh > ${work_dir}/complexes/prep.sh
#prepare free energy calculation script
sed "s#<RECEPTOR>#${receptor}# ; s#<LIGPAIRS>#${ligpreppairs}# " < ${work_dir}/scripts/free_energies_template.sh > ${work_dir}/complexes/free_energies.sh
