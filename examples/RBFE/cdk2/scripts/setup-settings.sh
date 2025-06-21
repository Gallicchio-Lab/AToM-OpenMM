#basename for jobs
basename=cdk2

#basename of the receptor pdb file (processed for amber using pdb4amber)
receptor=cdk2

#basenames of the ligand pairs (.sdf format expected)

ligands=(  "H1Q           H1R"
           "H1Q           H1S"
           "H1Q           OI9"
           "H1Q           OIU"
           "H1Q           OIY"
           "H1R           OIU"
           "OI9           H1S"            			
	   "OI9           OIY"
        )                                                                

ref_atoms=("14,21,18 14,21,18"   
           "14,21,18 15,22,19"
           "14,21,18 15,22,19"
           "14,21,18 16,23,20"
           "14,21,18 14,21,18"  
           "14,21,18 16,23,20"
           "15,22,19 15,22,19"
           "15,22,19 14,21,18")


#displacement vector
displacement=("-20.0" "0.0" "-20.0")

#residue ids of the receptor that define the center of the binding site 
#IMP::::::::::: same as what Meastro shows, setup-settings will do the math
vsite_rcpt_residues=( 12 14 16 22 84 87 88 134 146 147 )
