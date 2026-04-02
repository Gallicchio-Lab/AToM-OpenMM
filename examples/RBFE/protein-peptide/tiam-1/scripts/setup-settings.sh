#workflow mode:
#  ligand_rbfe       - existing ligand/SDF workflow
#  protein_mutation  - pair-specific PDB workflow using mutation residue backbone atoms
workflow_mode=protein_mutation

#basename of the receptor pdb file
receptor=tiam1

#reference ligand and reference alignment atoms
#used only in ligand_rbfe mode
ref_ligand=ligandBaseName
ref_ligand_alignment_atoms="14,21,18"

#transformations
#ligand_rbfe mode:       "LIG1 LIG2"
# not used in this example (workflow_mode=protein_mutation)
# see examples/RBFE/cdk2 for a complete ligand_rbfe example
ligand_transformations=()

#protein_mutation mode:  "MUT1 MUT2 RESID"
protein_mutation_transformations=(
    "sdc1E5Q sdc1wt 5"
    "sdc1wt sdc1E5Q 5"
    "sdc1Q3E sdc1wt 3"
    "sdc1wt sdc1Q3E 3"
)

case "${workflow_mode}" in
    ligand_rbfe)
        ligands=("${ligand_transformations[@]}")
        ;;
    protein_mutation)
        ligands=("${protein_mutation_transformations[@]}")
        ;;
    *)
        echo "Unsupported workflow_mode '${workflow_mode}' in setup-settings.sh" >&2
        return 1 2>/dev/null || exit 1
        ;;
esac
