import os
import argparse
import pickle
import yaml
import numpy as np
from pathlib import Path
from atom_openmm.utils.AtomUtils import get_alignment_atoms

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
    alignments = get_alignment_atoms(ref_lig_file, ref_lig_alignment_atoms, lig_files)

    #writes the alignments to a yaml file
    with open(args['alignmentsYAMLoutFile'], 'w') as file:
        yaml.dump(alignments, file, default_flow_style=None, width=1000000, sort_keys=False)
