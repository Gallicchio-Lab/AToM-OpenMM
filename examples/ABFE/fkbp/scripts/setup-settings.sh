
#basename for jobs
basename=fkbp

#basename of the receptor pdb file (processed for amber using pdb4amber)
receptor=fkbp

#basenames of the ligands (.sdf format expected)
ligands=(thi prp dss dmso dapp dap but)

#displacement vector
displacement=("22.0" "22.0" "22.0")

#residue ids of the receptor that define the center of the binding site
vsite_rcpt_residues=(26 36 37 42 46 48 54 55 56 82)
