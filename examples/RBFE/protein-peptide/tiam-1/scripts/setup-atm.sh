#!/bin/bash

work_dir=$(pwd)
scripts_dir=${work_dir}/scripts

#basename of the receptor pdb file
receptor=tiam1

#generate ligand indexes
(
    cd ${work_dir}/ligands || exit 1
    python ${scripts_dir}/make_refvar_atoms.py --mutations "sdc1wt sdc1A0F A8F, sdc1A0F sdc1wt F8A, sdc1wt sdc1A0M A8M, sdc1A0M sdc1wt M8A, sdc1wt sdc1A0V A8V, sdc1A0V sdc1wt A8V, sdc1wt sdc1A0V V8A, sdc1A0F sdc1A0M F8M, sdc1A0M sdc1A0F M8F, sdc1A0F sdc1A0V F8V, sdc1A0V sdc1A0F V8F, sdc1A0M sdc1A0V M8V, sdc1A0V sdc1A0M V8M" > ${scripts_dir}/ligand_indexes.sh || exit 1
)

. ${scripts_dir}/ligand_indexes.sh 

#displacement vector
displacement=("0.0" "40.0" "0.0")

#residue ids of the receptor that define the center of the binding site 
vsite_rcpt_residues=( 850 856 857 858 859 860 861 862 912 915 916 920 )


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

    #skip if already prepared
    if [ -d ${work_dir}/complexes/${jobname} ] ; then
	continue
    fi
    
    mkdir -p ${work_dir}/complexes/${jobname} || exit 1
    cd ${work_dir}/complexes/${jobname}    || exit 1

    #files of receptor and ligands
    rcptpdb=${work_dir}/receptor/${receptor}.pdb
    lig1pdb=${work_dir}/ligands/${lig1}.pdb
    lig2pdb=${work_dir}/ligands/${lig2}.pdb
    
    #creates simulation system
    displs=""
    for d in ${displacement[@]}; do
	displs="$displs $d"
    done
    echo "Displacement vector: $displs"
    make_atm_system_from_rcpt_lig --receptorinFile "${rcptpdb}" --LIG1inFile "${lig1pdb}" --LIG2inFile "${lig2pdb}" --displacement "${displs}" --systemXMLoutFile "${jobname}_sys.xml" --systemPDBoutFile "${jobname}.pdb" --hmass 1.5 --forcefieldJSONCachefile "${work_dir}/ligands/ffdb.json" || exit 1

    #ligand alignment reference atoms
    refpair=${ref_atoms[$l]}
    read -ra ref <<< $refpair
    ref_atoms1a="${ref[0]}"
    ref_atoms2a="${ref[1]}"

    #vsite residues
    vres=""
    for r in ${vsite_rcpt_residues[@]}; do
	r1=`expr $r - 0`
	vres="$vres $r1"
    done

    
    #ligand variable regions
    varregions=${variable_region[$ligpair]}
    read -ra var <<< $varregions
    var_lig1="${var[0]}"
    var_lig2="${var[1]}"

    lig1chainname=L
    lig2chainname=M
    
    cp ${work_dir}/scripts/vmd_template.in ${work_dir}/complexes/${jobname}/ || exit 1

    echo "python $scripts_dir/create_cntlfile_from_template_rbfe.py --systemPDBFile "${jobname}.pdb" --templatein ${work_dir}/scripts/asyncre_template.cntl --jobname "${jobname}" --displacement "${displs}" --lig1chainname "${lig1chainname}" --lig1refatoms "${ref_atoms1a}" --lig2chainname "${lig2chainname}" --lig2refatoms "${ref_atoms2a}" --lig1varatoms "${var_lig2}"--lig1varatoms "${var_lig2}"  --vsiteResidues "${vres}" --cntlfileout ${jobname}_asyncre.cntl "
    
    python $scripts_dir/create_cntlfile_from_template_rbfe.py --systemPDBFile "${jobname}.pdb" --templatein ${work_dir}/scripts/asyncre_template.cntl --jobname "${jobname}" --displacement "${displs}" --lig1chainname "${lig1chainname}" --lig1refatoms "${ref_atoms1a}" --lig2chainname "${lig2chainname}" --lig2refatoms "${ref_atoms2a}" --lig1varatoms "${var_lig1}" --lig2varatoms "${var_lig2}" --vsiteResidues "${vres}" --cntlfileout ${jobname}_asyncre.cntl || exit 1

    sed "s#<JOBNAME>#${jobname}#g;s#<CONDAENV>#${CONDA_DEFAULT_ENV}#" < ${work_dir}/scripts/run_template.sh > ${work_dir}/complexes/${jobname}/run.sh || exit 1
    sed "s#<JOBNAME>#${jobname}#g;s#<CONDAENV>#${CONDA_DEFAULT_ENV}#" < ${work_dir}/scripts/run4_template.sh > ${work_dir}/complexes/${jobname}/run4.sh || exit 1

    cp ${work_dir}/scripts/analyze.sh ${work_dir}/scripts/uwham_analysis.R ${work_dir}/complexes/${jobname}/ || exit 1
    
done

#prepare free energy calculation script
sed "s#<RECEPTOR>#${receptor}# ; s#<LIGPAIRS>#${ligpreppairs}# " < ${work_dir}/scripts/free_energies_template.sh > ${work_dir}/complexes/free_energies.sh
