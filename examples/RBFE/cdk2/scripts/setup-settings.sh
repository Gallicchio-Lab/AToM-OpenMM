#basename for jobs
basename=cdk2

#path to ASyncRE
asyncre_dir=~/AToM-OpenMM

#basename of the receptor pdb file (processed for amber using pdb4amber)
receptor=cdk2
#number of atoms of the protein receptor not including water

#basenames of the ligand pairs (.mol2 format expected)
#with their reference alignment atoms

#1h1q 14,21,18

#1h1r 14,21,18

#1h1s 15,22,19

#1oiu 16,23,20

#1oiy 14,21,18

#1oi9 15,22,19


ligands=(  "1H1Q           1H1R"
           "1H1Q           1H1S"
           "1H1Q           1OI9"
           "1H1Q           1OIU"
           "1H1Q           1OIY"
           "1H1R           1OIU"
           "1OI9           1H1S"            			
	   "1OI9           1OIY"
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

