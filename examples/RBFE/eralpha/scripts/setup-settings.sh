
#basename for jobs
basename=eralpha

#basename of the receptor pdb file (processed for amber using pdb4amber)
receptor=eralpha

#basenames of the ligand pairs (.mol2 format expected)
#with their reference alignment atoms

  ligands=("2d      2e"      "2d      3a"      "2d      3b"      "2e      3a"      "2e      3b"      "3a      3b")
ref_atoms=("7,11,23 12,7,21" "7,11,23 15,10,5" "7,11,23 13,8,22" "12,7,21 15,10,5" "12,7,21 13,8,22" "15,10,5 13,8,22")

#displacement vector
displacement=("22.0" "22.0" "22.0")

#residue ids of the receptor that define the center of the binding site
vsite_rcpt_residues=(36 39 40 42 43 77 80 81 84 95 97 111 113 114 117 118 214 215 217 218)

