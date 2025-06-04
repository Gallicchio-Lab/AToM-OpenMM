
#basename for jobs
basename=TYK2-JAK2

#path to AToM
AToM_dir=$HOME/AToM-OpenMM/AToM-OpenMM

#Maestro-prepared receptor pdb files
receptor1=TYK2
receptor2=JAK2

#ligand pairs (.sdf format expected)
#with their reference alignment atoms
ligands=("c46	c23"
		"c23 	c46")

ref_atoms=("5,12,2 		15,1,7"
			"15,1,7   	5,12,2")

#displacement vector
displacement=("0.0" "70.0" "0.0")

#residue ids of the receptor that define the center of the binding site
vsite_rcpt_residues=(892)

