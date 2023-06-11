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
    python $AToM_dir/make_atm_system_from_rcpt_lig.py --receptorinFile "${rcptpdb}" --LIG1SDFinFile "${ligsdf}" --displacement "${displs}" --systemXMLoutFile "${jobname}_sys.xml" --systemPDBoutFile "${jobname}.pdb" --forcefieldJSONCachefile 'ligands/ffdb.json' || exit 1

    #convert list of resids of receptor into Calpha atom indexes
    #these atoms are restrained.
    #Indexes are shifted by 1 as they are expected to start from 0 in the OpenMM scripts
    calpha_atoms=$( awk '/^ATOM/ && $3 ~ /^CA$/ {print $2}' ${jobname}.pdb ) || exit 1
    unset restr_atoms
    for i in $calpha_atoms; do
	i1=$( expr $i - 1 )
	if [ -z "$restr_atoms" ] ; then
	    restr_atoms=$i1
	else
	    restr_atoms="${restr_atoms}, $i1"
	fi
    done
    echo "Indexes of the Calpha atoms:"
    echo "$restr_atoms"

    #get the list of receptor atoms that define the centroid of the binding site
    calphascm=()
    n=0
    for res in ${vsite_rcpt_residues[@]}; do 
	i=$( awk -v resid=$res '/^ATOM/ && $6 == resid && $3 ~ /^CA$/ {print $2}' ${jobname}.pdb ) || exit 1
	calphascm[$n]=$i
	n=$(expr $n + 1)
    done
    unset vsite_rcpt_atoms
    for i in ${calphascm[@]}; do
	i1=$( expr $i - 1 )
	if [ -z "$vsite_rcpt_atoms" ] ; then
	    vsite_rcpt_atoms=$i1
	else
	    vsite_rcpt_atoms="${vsite_rcpt_atoms}, $i1"
	fi
    done
    echo "Vsite receptor atoms:"
    echo "$vsite_rcpt_atoms"
 
    #search for the first ligand to find the size of the receptor
    tlig=${lig:0:2}
    namereslig1=${tlig^^}
    l1=$( awk "\$4 ~ /${namereslig1}/{print \$2}" < ${jobname}.pdb  | head -1 ) || exit 1
    lig1start=$(expr $l1 - 1 )
    num_atoms_rcpt=$(expr $l1 - 1 )
    #retrieve number of atoms of the ligands from their mol2 files
    num_atoms_lig1=$( awk 'NR==4{print $1}' ${work_dir}/ligands/${lig}.sdf ) || exit 1
    
    #list of ligand atoms (starting from zero) assuming they are listed soon after the receptor
    lig1end=$( expr $lig1start + $num_atoms_lig1 - 1 )
    lig1_atoms=""
    for i in $( seq $lig1start $lig1end ); do
	if [ -z "$lig1_atoms" ] ; then
	    lig1_atoms=$i
	else
	    lig1_atoms="${lig1_atoms}, $i"
	fi
    done
    
    echo "atoms of $lig are $lig1_atoms"
    
    displs=""
    count=0
    for d in ${displacement[@]}; do
	displs="$displs $d"
	count=$(expr $count + 1 )
    done
    echo "Displacement vector: $displs"

    #builds asyncre scripts
    replstring="s#<JOBNAME>#${jobname}# ; s#<DISPLX>#${displacement[0]}# ; s#<DISPLY>#${displacement[1]}# ; s#<DISPLZ>#${displacement[2]}# ; s#<VSITERECEPTORATOMS>#${vsite_rcpt_atoms}# ; s#<RESTRAINEDATOMS>#${restr_atoms}# ; s#<LIGATOMS>#${lig1_atoms}#" 
    sed "${replstring}" < ${work_dir}/scripts/asyncre_template.cntl > ${jobname}_asyncre.cntl || exit 1

    #copy slurm files, etc
    sed "s#<JOBNAME>#${jobname}#;s#<ASYNCRE_DIR>#${AToM_dir}#" < ${work_dir}/scripts/run_template.sh > ${work_dir}/complexes/${jobname}/run.sh
    cp ${work_dir}/scripts/analyze.sh ${work_dir}/scripts/uwham_analysis.R ${work_dir}/complexes/${jobname}/
    
done

#prepare prep script
ligs=${ligands[@]}
sed "s#<RECEPTOR>#${receptor}#;s#<LIGS>#${ligs}#  ; s#<ASYNCRE_DIR>#${AToM_dir}#g "< ${work_dir}/scripts/prep_template.sh > ${work_dir}/complexes/prep.sh

#prepare free energy calculation script
sed "s#<RECEPTOR>#${receptor}#;s#<LIGS>#${ligs}# " < ${work_dir}/scripts/free_energies_template.sh > ${work_dir}/complexes/free_energies.sh
