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

    #assign GAFF2 parameters to ligands
    cd ${work_dir}/ligands || exit 1
    for lig in $lig1 $lig2 ; do
	if [ ! -f ${lig}-p.mol2 ] ; then
	    charge=$( awk 'BEGIN{charge = 0} ; NF == 9 {charge += $9} ; END {if (charge >= 0) {print int(charge + 0.5)} else {print int(charge - 0.5)}}' < ${lig}.mol2 ) || exit 1
	    echo "antechamber -pl 15 -fi mol2 -fo mol2 -i ${lig}.mol2 -o ${lig}-p.mol2 -c bcc -at gaff2 -nc ${charge} "
	    antechamber -pl 15 -fi mol2 -fo mol2 -i ${lig}.mol2 -o ${lig}-p.mol2 -c bcc -at gaff2 -nc ${charge} || exit 1
	    echo "parmchk2 -i ${lig}-p.mol2 -o ${lig}-p.frcmod -f mol2"
	    parmchk2 -i ${lig}-p.mol2 -o ${lig}-p.frcmod -f mol2 || exit 1
	fi
    done

    #creates system in complexes folder
    jobname=${receptor}-${lig1}-${lig2}
    mkdir -p ${work_dir}/complexes/${jobname} || exit 1
    cd ${work_dir}/complexes/${jobname}    || exit 1
    
    displs=""
    for d in ${displacement[@]}; do
	displs="$displs $d"
    done
    echo "Displacement vector: $displs"

    (
cat <<EOF
  source leaprc.protein.ff14SB
  source leaprc.gaff2
  source leaprc.water.tip3p
  RCPT = loadpdb "${work_dir}/receptor/${receptor}.pdb"
  LIG1 = loadmol2 "${work_dir}/ligands/${lig1}-p.mol2"
  loadamberparams "${work_dir}/ligands/${lig1}-p.frcmod"
  LIG2 = loadmol2 "${work_dir}/ligands/${lig2}-p.mol2"
  loadamberparams "${work_dir}/ligands/${lig2}-p.frcmod"
  translate LIG2 { ${displs}  }
  MOL = combine {RCPT LIG1 LIG2}
  addions2 MOL Na+ 0
  addions2 MOL Cl- 0
  solvateBox MOL TIP3PBOX 10.0
  saveamberparm MOL ${jobname}.prmtop ${jobname}.inpcrd
  quit
EOF
) > tleap.cmd

    tleap -f tleap.cmd || exit 1
    
    #convert prmtop/inpcrd to OpenMM's System + PDB topology
    make_atm_system_from_amber --AmberPrmtopinFile ${jobname}.prmtop --AmberInpcrdinFile ${jobname}.inpcrd --systemXMLoutFile ${jobname}_sys.xml --systemPDBoutFile ${jobname}.pdb

    #residue ligand names
    lig1resname="UNK"
    lig2resname="UNK"
    
    #ligand alignment reference atoms
    refpair=${ref_atoms[$l]}
    read -ra ref <<< $refpair
    ref_atoms1a="${ref[0]}"
    ref_atoms2a="${ref[1]}"

    #vsite residues
    vres=${vsite_rcpt_residues[@]}

    echo "python $scripts_dir/create_cntlfile_from_template_rbfe.py --systemPDBFile "${jobname}.pdb" --templatein ${work_dir}/scripts/asyncre_template.cntl --jobname "${jobname}" --displacement "${displs}" --lig1resname "${lig1resname}" --lig1refatoms "${ref_atoms1a}" --lig2resname "${lig2resname}" --lig2refatoms "${ref_atoms2a}"  --vsiteResidues "${vres}" --cntlfileout ${jobname}_asyncre.cntl "
    
    python $scripts_dir/create_cntlfile_from_template_rbfe.py --systemPDBFile "${jobname}.pdb" --templatein ${work_dir}/scripts/asyncre_template.cntl --jobname "${jobname}" --displacement "${displs}" --lig1resname "${lig1resname}" --lig1refatoms "${ref_atoms1a}" --lig2resname "${lig2resname}" --lig2refatoms "${ref_atoms2a}"  --vsiteResidues "${vres}" --cntlfileout ${jobname}_asyncre.cntl || exit 1
    
    
    #copy slurm files, etc
    sed "s#<JOBNAME>#${jobname}#" < ${work_dir}/scripts/run_template.sh > ${work_dir}/complexes/${jobname}/run.sh

    cp ${work_dir}/scripts/analyze.sh ${work_dir}/scripts/uwham_analysis.R ${work_dir}/complexes/${jobname}/
    
done

#prepare prep script
sed "s#<RECEPTOR>#${receptor}# ; s#<LIGPAIRS>#${ligpreppairs}#g " < ${work_dir}/scripts/prep_template.sh > ${work_dir}/complexes/prep.sh
#prepare free energy calculation script
sed "s#<RECEPTOR>#${receptor}# ; s#<LIGPAIRS>#${ligpreppairs}# " < ${work_dir}/scripts/free_energies_template.sh > ${work_dir}/complexes/free_energies.sh
