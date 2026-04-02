import os
import argparse
import pickle
import yaml
import numpy as np
from scipy.spatial import distance_matrix
from scipy.spatial.distance import cdist
from pathlib import Path
from rdkit import Chem

def _get_solute_coords(solute_fpath: Path):  # pdb or sdf file format
    """
    Return N*3 array for solute coordinates
    """
    assert solute_fpath.suffix in [".sdf", ".pdb"], f"{solute_fpath} is not supported."
    if solute_fpath.suffix == ".sdf":
        mol = Chem.SDMolSupplier(str(solute_fpath), removeHs=False)[0]
    else:
        mol = Chem.rdmolfiles.MolFromPDBFile(str(solute_fpath), removeHs=False)

    conf = mol.GetConformer()
    N_atoms = mol.GetNumAtoms()
    coords = np.zeros((N_atoms, 3))
    for row in range(N_atoms):
        coords[row] = np.array(list(conf.GetAtomPosition(row)))
    return coords

def get_alignment(ref_lig_file, ref_lig_alignment_atoms, lig_files):
    """
    Return Dict with key = 'ligand_name`,
        value= dict {'aligned_atom_ids': [id-1,id_2,id_3], 'N_atoms': N}
    """
    lig_paths = [Path(f) for f in lig_files]
    ligand_names = [path.stem for path in lig_paths]

    result = {}

    #dictionary of ligand coordinates
    ligands_coords = { p.stem: _get_solute_coords(p) for p in lig_paths }
    
    ref_lig_path = Path(ref_lig_file)
    ref_lig_coords = _get_solute_coords(ref_lig_path)
    ref_lig_N_atoms = ref_lig_coords.shape[0]
    ref_lig_name = ref_lig_path.stem

    result.update(
        {
            ref_lig_name: {
                "align_atom_ids": ref_lig_alignment_atoms,
                "N_atoms": ref_lig_N_atoms,
            }
        }
    )

    for lig_name, lig_coords in ligands_coords.items():
        dist_matrix = cdist(ref_lig_coords, lig_coords)
        nearest_ids = np.argmin(dist_matrix, axis=1)
        # the ref_align_idx starts from 1, be careful about the index
        lig_aligment_ids = [nearest_ids[i - 1] + 1 for i in ref_lig_alignment_atoms]
        result.update(
            {
                lig_name: {
                    "align_atom_ids": [int(i) for i in lig_aligment_ids],
                    "N_atoms": lig_coords.shape[0],
                }
            }
        )
    return result

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    #required arguments
    parser.add_argument('--refligFile',        help='The structure file of the reference ligand usually in SDF format', required=True)
    parser.add_argument('--LIG1refatoms',      help='The alignment atoms of the reference ligand, --LIG1refatoms "14,12,11" for example, starting with 1', type=str, required=True)
    parser.add_argument('--alignmentsYAMLoutFile', help='The yaml-formatted file where to save the alignments.', type=str, required=True)
    #optional arguments
    parser.add_argument('--ligFilesGlob',      help='A glob pattern to find all the other ligands', default='*.sdf')

    #dictionary of arguments
    args = vars(parser.parse_args())

    #list of ligand files
    p = Path('.')
    # .glob() returns a generator (efficient for thousands of files)
    # use list() if you need to sort them or access them by index
    ligfiles_path = list(p.glob(args['ligFilesGlob']))
    #to access the file names:
    lig_files = [f.name for f in ligfiles_path]
    
    ref_lig_file = args['refligFile']
    ref_lig_alignment_atoms = [int(i) for i in args['LIG1refatoms'].split(',')]
    alignments = get_alignment(ref_lig_file, ref_lig_alignment_atoms, lig_files)

    #writes the alignments to a yaml file
    with open(args['alignmentsYAMLoutFile'], 'w') as file:
        yaml.dump(alignments, file, default_flow_style=None, width=1000000, sort_keys=False)
