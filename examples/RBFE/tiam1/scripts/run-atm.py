import os
import argparse
import yaml
import numpy as np
from pathlib import Path
import pandas as pd
import openmm as mm
from openmm import *
from openmm.app import *
from openmm.unit import *
#from multiprocessing import freeze_support
from atom_openmm.make_atm_system_from_rcpt_lig import make_system
from atom_openmm.rbfe_structprep import rbfe_structprep
from atom_openmm.rbfe_production import rbfe_production
from atom_openmm.utils.AtomUtils import AtomUtils, get_selected_principal_groups
from atom_openmm.uwham import calculate_uwham_from_rundir
#from openmm_run import openmm_run

PROTEIN_BACKBONE_ATOM_NAMES = {
    "N", "NT", "CA", "CAY", "CAT", "C", "CY", "O", "OY", "OXT",
    "H", "H2", "H3", "HY1", "HY2", "HY3", "HA",
    "DN", "DCA", "DC", "DO",
    "LPOA", "LPOB", "DOT1", "DOT2", "LPT1", "LPT2", "LPT3", "LPT4",
}

def get_indexes_from_query(topology, query):
    indexes = [ atom.index for atom in topology.atoms() if eval(query, {"atom": atom}) ]
    indexes.sort()
    return indexes

#get the indexes of the atoms of a residue. Optionally, filter them by a query
def get_indexes_from_residue(residue, query = 'True'):
    indexes = [
        atom.index
        for atom in residue.atoms()
        if eval(query, {"atom": atom, "PROTEIN_BACKBONE_ATOM_NAMES": PROTEIN_BACKBONE_ATOM_NAMES})
    ]
    indexes.sort()
    return indexes

def get_indexes_from_residues(residues, query='True'):
    indexes = []
    for residue in residues:
        indexes.extend(
            atom.index
            for atom in residue.atoms()
            if eval(query, {"atom": atom, "PROTEIN_BACKBONE_ATOM_NAMES": PROTEIN_BACKBONE_ATOM_NAMES})
        )
    indexes.sort()
    return indexes

def get_residues_by_name(topology, residue_name):
    residues = []
    for chain in topology.chains():
        for residue in chain.residues():
            if residue.name == residue_name:
                residues.append(residue)
    return residues

def get_residues_by_chain_id(topology, chain_id):
    residues = []
    for chain in topology.chains():
        if chain.id == chain_id:
            residues.extend(list(chain.residues()))
    return residues

def get_residue_by_id(residues, residue_id):
    matches = [residue for residue in residues if residue.id == str(residue_id)]
    if not matches:
        return None
    if len(matches) > 1:
        raise ValueError(
            f"Found multiple residues with id {residue_id} in the selected alchemical partner."
        )
    return matches[0]

#calculates the center of a set of atoms
def cm_from_indexes(topology, positions, indexes):
    cm = Vec3(0,0,0)*nanometer
    n = 0
    for atom in topology.atoms():
        if atom.index in indexes:
            cm += positions[atom.index]
            n += 1
    cm = cm/float(n)
    return cm

def _get_solute_coords(solute_fpath: Path):  # pdb or sdf file format
    """
    Return N*3 array for solute coordinates
    """
    assert solute_fpath.suffix in [".sdf", ".pdb"], f"{solute_fpath} is not supported."
    if solute_fpath.suffix == ".sdf":
        from rdkit import Chem
        mol = Chem.SDMolSupplier(str(solute_fpath), removeHs=False)[0]
    else:
        coords = []
        with open(solute_fpath, "r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    coords.append(
                        [
                            float(line[30:38]),
                            float(line[38:46]),
                            float(line[46:54]),
                        ]
                    )
        return np.asarray(coords, dtype=float)

    conf = mol.GetConformer()
    coords = np.zeros((mol.GetNumAtoms(), 3))
    for row in range(mol.GetNumAtoms()):
        coords[row] = np.array(list(conf.GetAtomPosition(row)))
    return coords

#adapted Eric Chen's atm implementation https://github.com/EricChen521/atm
def calc_displ_vec(receptor_file, ligand2_file, options):
    """
    Return the optimal displacement vector (x,y,z) for the second ligand.

    Step 1: Find the smallest area center(point_1),
    Step 2: Move point_1 along the third direction with 9 A to point_2,
    Step 3: find the the smallest point in the third direction, point_3,
    Finally, the displacement_vec is obtained point_2 - point_3
    """
    protein_fpath = Path(receptor_file)
    ligand_dpath = Path(ligand2_file)

    protein_coords = _get_solute_coords(solute_fpath=protein_fpath)
    x_range = np.array([0, min(protein_coords[:, 0]), max(protein_coords[:, 0])])
    y_range = np.array([1, min(protein_coords[:, 1]), max(protein_coords[:, 1])])
    z_range = np.array([2, min(protein_coords[:, 2]), max(protein_coords[:, 2])])

    # print(f"system size: X {x_range}, Y {y_range}, Z {z_range}")

    small_area_center = np.array([0.0, 0.0, 0.0])

    axes = sorted([x_range, y_range, z_range], key=lambda v: v[2] - v[1])
    small_area_center[int(axes[0][0])] = np.mean(axes[0][1:])
    small_area_center[int(axes[1][0])] = np.mean(axes[1][1:])
    small_area_center[int(axes[2][0])] = axes[2][2] # max u coordinate
    u = int(axes[2][0])
    # point_1
    #print(f"small_area_center coordiante: {small_area_center}, with u axis: {u}")
    # point_2
    small_area_center[u] += 10.0 # max u coordinate + 10

    # find the smallest U in all ligands.
    ligands_coords =_get_solute_coords(ligand_dpath)

    sorted_u_index = ligands_coords[:, u].argsort()
    sorted_u_coords = ligands_coords[sorted_u_index]
    distant_lig_atom_coords = sorted_u_coords[0, :] # minimum u of all ligands
    #print(f"sorted ligand atom coordiates by {u}: {sorted_u_coords}")

    displ_vec = np.round((small_area_center - distant_lig_atom_coords), 2)
    print(f"Automatic displacement vector: {displ_vec}")
    return displ_vec

def rbfe_setup(receptor_file, lig1_file, lig2_file, ff_json_file, options):
    #prepare the system, unless already done
    basename = options['BASENAME']
    setup = {}
    setup['receptorfile'] = receptor_file
    setup['lig1file'] = lig1_file
    setup['lig2file'] = lig2_file
    setup['displacement'] = options['DISPLACEMENT']
    setup['xmloutfile'] = basename + "_sys.xml"
    setup['pdboutfile'] = basename + ".pdb"
    setup['ffcachefile'] = ff_json_file
    if 'HMASS' in options:
        setup['hmass'] = options['HMASS']
    if 'LIGAND_FORCE_FIELD' in options:
        setup['ligandforcefield'] = options['LIGAND_FORCE_FIELD']
    #add more later ...

    make_system(**setup)

def rbfe_prepare_args(options):
    basename = options['BASENAME']
    
    #system from pdb file
    pdb = PDBFile(basename + ".pdb")
    topology = pdb.topology
    positions = pdb.positions

    #LIGAND1_ATOMS (prefer residues named L1; fallback to chain L)
    lig1resname = 'L1'
    res1 = get_residues_by_name(topology, lig1resname)
    if not res1:
        res1 = get_residues_by_chain_id(topology, 'L')
    if not res1:
        raise ValueError("Could not identify alchemical partner 1. Expected residue name L1 or chain ID L in the prepared system.")
    ligand1_atoms = get_indexes_from_residues(res1)
    options['LIGAND1_ATOMS'] = ligand1_atoms

    #LIGAND2_ATOMS (prefer residues named L2; fallback to chain M)
    lig2resname = 'L2'
    res2 = get_residues_by_name(topology, lig2resname)
    if not res2:
        res2 = get_residues_by_chain_id(topology, 'M')
    if not res2:
        raise ValueError("Could not identify alchemical partner 2. Expected residue name L2 or chain ID M in the prepared system.")
    ligand2_atoms = get_indexes_from_residues(res2)
    options['LIGAND2_ATOMS'] = ligand2_atoms

    # In protein_mutation mode, only the mutation-site sidechains are alchemically variable.
    if options.get('WORKFLOW_MODE') == 'protein_mutation':
        mutation_residue = options.get('PAIR_ALIGNMENT_RESIDUE')
        if not mutation_residue:
            raise ValueError(
                "protein_mutation mode requires PAIR_ALIGNMENT_RESIDUE to define variable atoms."
            )

        res1_var = get_residue_by_id(res1, mutation_residue)
        if res1_var is None:
            raise ValueError(
                f"Could not find mutation residue {mutation_residue} in alchemical partner 1."
            )
        res2_var = get_residue_by_id(res2, mutation_residue)
        if res2_var is None:
            raise ValueError(
                f"Could not find mutation residue {mutation_residue} in alchemical partner 2."
            )

        options['LIGAND1_VAR_ATOMS'] = get_indexes_from_residue(
            res1_var, query='atom.name not in PROTEIN_BACKBONE_ATOM_NAMES'
        )
        options['LIGAND2_VAR_ATOMS'] = get_indexes_from_residue(
            res2_var, query='atom.name not in PROTEIN_BACKBONE_ATOM_NAMES'
        )
    else:
        #the whole ligands are variable (dual topology)
        options['LIGAND1_VAR_ATOMS'] = options['LIGAND1_ATOMS']
        options['LIGAND2_VAR_ATOMS'] = options['LIGAND2_ATOMS']

    #ligand anchor atoms
    ligand1_ref_atoms = options['ALIGN_LIGAND1_REF_ATOMS']
    ligand2_ref_atoms = options['ALIGN_LIGAND2_REF_ATOMS']
    lig1_anchor_atom = ligand1_atoms[ligand1_ref_atoms[0]]
    lig2_anchor_atom = ligand2_atoms[ligand2_ref_atoms[0]]
    options['LIGAND1_ATTACH_ATOM'] =  lig1_anchor_atom
    options['LIGAND2_ATTACH_ATOM'] =  lig2_anchor_atom

    #ligand CM atoms (same as anchor atoms)
    options['LIGAND1_CM_ATOMS'] = [ ligand1_atoms[ ligand1_ref_atoms[0] ] ]
    options['LIGAND2_CM_ATOMS'] = [ ligand2_atoms[ ligand2_ref_atoms[0] ] ]

    #position of CMs of the ligands
    lig1cm_pos = cm_from_indexes(topology, positions, options['LIGAND1_CM_ATOMS'] )
    lig2cm_pos = cm_from_indexes(topology, positions, options['LIGAND2_CM_ATOMS'] )

    #displacement (distance between the ligand CMs)
    displ = (lig2cm_pos - lig1cm_pos).value_in_unit(angstrom)
    options['DISPLACEMENT'] = [ displ.x , displ.y, displ.z ] 

    #CM atoms of the receptor, all C-alpha atoms
    rcpt_chain_name = options.get('RCPT_CHAIN_NAME', 'A')
    rcpt_frame_query = f'atom.residue.chain.id == "{rcpt_chain_name}" and atom.name == "CA"'
    rcpt_frame_indexes = get_indexes_from_query(topology, rcpt_frame_query)

    #internal receptor frame
    rcpt_frame = get_selected_principal_groups(topology, positions, rcpt_frame_indexes)

    #receptor CM atoms (origin of the internal frame)
    options['RCPT_CM_ATOMS'] = rcpt_frame['origin']['indices']
    options['RCPT_FRAME_ATOMS_O'] = options['RCPT_CM_ATOMS']
    options['RCPT_FRAME_ATOMS_Z'] = rcpt_frame['z_axis']['indices']
    options['RCPT_FRAME_ATOMS_Y'] = rcpt_frame['y_axis']['indices']

    #offset is the distance from the CM of the first ligand from the CM of the receptor
    rcpt_cm_pos = Vec3(float(rcpt_frame['origin']['com'][0]),
                       float(rcpt_frame['origin']['com'][1]),
                       float(rcpt_frame['origin']['com'][2])) * nanometer
    offset = (lig1cm_pos - rcpt_cm_pos).value_in_unit(angstrom)
    options['LIGOFFSET'] = [ offset.x, offset.y, offset.z ]

    #restrained atoms (same as RCPT_CM_ATOMS)
    options['POS_RESTRAINED_ATOMS'] = options['RCPT_CM_ATOMS']

    #receptor-ligand exclusion potential
    options['EXCLUSION_POT_MOL1_INDEXES'] = get_indexes_from_query(
        topology, f'(atom.residue.chain.id == "{rcpt_chain_name}") and (atom.element.atomic_number != 1)')
    options['EXCLUSION_POT_MOL2_INDEXES'] = get_indexes_from_residues(
        res2, query='(atom.element.atomic_number != 1)'
    )

    #working directory
    options['WORKDIR'] = os.environ.get('PWD') 

def rbfe_production_data(options):
    rundir = options['WORKDIR']
    jobname = options['BASENAME']
    nstates = len(options['LAMBDAS'])
    # Create list of data files
    datafiles = [
        os.path.join(rundir, f"r{i}", f"{jobname}.out") for i in range(nstates)
    ]
    # find the smallest number of samples for each state
    line_counts = {}
    for file_path in datafiles:
        try:
            count = 0
            with pd.read_csv(file_path, chunksize=10000, header=None, sep=None, engine='python') as reader:
                for chunk in reader:
                    count += len(chunk)
            line_counts[file_path] = count
        except:
            pass
    # Find the the minimum count
    if line_counts:
        min_file = min(line_counts, key=line_counts.get)
        min_value = line_counts[min_file]
        return min_value
    else:
        return None

def create_vmd_infile(options):
    template_path = 'vmd_template.in'
    output_path = 'vmd.in'
    replacements = {
        "<LIG1ATTACHATOM>": options['LIGAND1_ATTACH_ATOM'],
        "<LIG2ATTACHATOM>": options['LIGAND2_ATTACH_ATOM']
    }
    with open(template_path, 'r') as f:
        content = f.read()
    # Loop through the dictionary and apply each replacement
    for key, value in replacements.items():
        content = content.replace(key, str(value))
    with open(output_path, 'w') as f:
        f.write(content)

def run_atm(options,
            jobname,
            receptor_file,
            lig1_file,
            lig2_file,
            workflow_mode='ligand_rbfe',
            alignments = None,
            pair_alignments = None,
            lig1_ref_atoms = None,
            lig2_ref_atoms = None,
            ff_json_file = 'ff.json'):

    #basename
    options['BASENAME'] = jobname
    options['WORKFLOW_MODE'] = workflow_mode

    #alignment atoms
    if workflow_mode == 'protein_mutation':
        assert(pair_alignments is not None)
        assert(jobname in pair_alignments)
        pair_alignment = pair_alignments[jobname]
        lig1_ref_atoms = pair_alignment['lig1_align_atom_ids']
        lig2_ref_atoms = pair_alignment['lig2_align_atom_ids']
        options['PAIR_ALIGNMENT_RESIDUE'] = pair_alignment.get('residue')
        options['PAIR_ALIGNMENT_ATOM_ORDER'] = pair_alignment.get('atom_order', ['CA', 'N', 'C'])
    elif alignments is not None:
        lig1_name = Path(lig1_file).stem
        lig2_name = Path(lig2_file).stem
        assert(lig1_name in alignments)
        assert(lig2_name in alignments)
        lig1_ref_atoms = alignments[lig1_name]['align_atom_ids']
        lig2_ref_atoms = alignments[lig2_name]['align_atom_ids']
    assert(lig1_ref_atoms is not None and lig2_ref_atoms is not None)
    options['ALIGN_LIGAND1_REF_ATOMS'] = [i-1 for i in lig1_ref_atoms]
    options['ALIGN_LIGAND2_REF_ATOMS'] = [i-1 for i in lig2_ref_atoms]

    #figures out an optimal displacement if one is not provided
    if not options.get('DISPLACEMENT'):
        displ = calc_displ_vec(receptor_file, lig2_file, options)
        options['DISPLACEMENT'] = list(displ)
    assert(options['DISPLACEMENT'])

    #sets up the system
    if not Path(jobname + '.pdb').exists():
        rbfe_setup(receptor_file, lig1_file, lig2_file, ff_json_file, options)

    #calculate atom indexes, distances, internal frame
    rbfe_prepare_args(options)

    #writes options to a yaml file for later executions
    if not Path(options['BASENAME'] + ".yaml").exists():
        with open(options['BASENAME'] + '.yaml', 'w') as file:
            yaml.dump(options, file, default_flow_style=None, width=1000000, sort_keys=False)

    #prepares vmd.in file if the template is present
    if Path('vmd_template.in').exists() and not Path('vmd.in').exists():
        create_vmd_infile(options)
        
    #minimization, equilibration, unless already done
    if not Path(options['BASENAME'] + "_0.xml").exists():
        rbfe_structprep(config_file = None, options = options)

    #assess existing production if any
    nsamples_per_replica = rbfe_production_data(options)
        
    #production
    if ( not options.get('MAX_SAMPLES', None) or 
         not nsamples_per_replica or
         nsamples_per_replica < options['MAX_SAMPLES'] ):
        rbfe_production(config_file = None, options = options)

    #assess final production if any
    nsamples_per_replica = rbfe_production_data(options)

    #free eneergy analysis
    if nsamples_per_replica:
        discard = int(nsamples_per_replica/3)
        ddG, ddG_std, uwham_data = calculate_uwham_from_rundir(
            options['WORKDIR'],
            options['BASENAME'],
            mintimeid=discard,
        )
        dgbind1 = uwham_data['dg_leg1']
        dgbind2 = uwham_data['dg_leg2']
        nsamples = uwham_data['nsamples']
        print(f'Relative binding free energy estimate after {nsamples+discard} perturbation energy samples per replica discarding the first {discard} samples:')
        print('DDGb=', ddG, '+/-', ddG_std, 'kcal/mol')
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    #required arguments
    parser.add_argument('--workflowMode',            help='Workflow mode: ligand_rbfe or protein_mutation', default='ligand_rbfe')
    parser.add_argument('--optionsYAMLinFile',        help='A YAML file with the default options (alchemical schedule, etc.)', required=True)
    parser.add_argument('--jobBasename',              help='The basename of the job', required=True)
    parser.add_argument('--receptorinFile',           help='The structure file of the receptor usually in PDB format', required=True)
    parser.add_argument('--LIG1inFile',               help='The structure file of the first alchemical partner, usually in SDF or PDB format', required=True)
    parser.add_argument('--LIG2inFile',               help='The structure file of the second alchemical partner, usually in SDF or PDB format', required=True)

    #optional arguments.
    #the alignment atoms can be specified explicitly or a pickle file with the alignment dictionary must be given
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument( '--alignmentsYAMLinFile',     help='The yaml file with the alignment dictionary', type=str)
    group.add_argument( '--pairAlignmentsYAMLinFile', help='The yaml file with the pair-specific alignment dictionary', type=str)
    group.add_argument( '--LIG1refatoms',             help='The alignment atoms of the first ligand, --LIG1refatoms "14,12,11" for example, starting with 1', type=str)
    parser.add_argument('--LIG2refatoms',             help='The alignment atoms of the second ligand, --LIG2refatoms "14,12,12" for example, starting with 1', type=str)

    #optional arguments
    parser.add_argument('--forcefieldJSONCachefile', help='file to store force field information', default='ff.json')
    
    #dictionary of arguments
    args = vars(parser.parse_args())

    jobname = args['jobBasename']
    lig1_file = args['LIG1inFile']
    lig2_file = args['LIG2inFile']
    receptor_file = args['receptorinFile']
    ff_json_file = args['forcefieldJSONCachefile']
    workflow_mode = args['workflowMode']
    
    #load default settings
    options = None
    with open(args['optionsYAMLinFile'], 'r') as file:
        options = yaml.safe_load(file)
    assert(options)

    alignments = None
    pair_alignments = None
    lig1_ref_atoms = None
    lig2_ref_atoms = None
    if workflow_mode not in ('ligand_rbfe', 'protein_mutation'):
        parser.error("--workflowMode must be one of: ligand_rbfe, protein_mutation")

    if workflow_mode == 'protein_mutation':
        if args['pairAlignmentsYAMLinFile'] is None:
            parser.error("--workflowMode protein_mutation requires --pairAlignmentsYAMLinFile.")
        with open(args['pairAlignmentsYAMLinFile'], 'r') as f:
            pair_alignments = yaml.safe_load(f)
    elif args['pairAlignmentsYAMLinFile'] is not None:
        parser.error("--pairAlignmentsYAMLinFile is only supported with --workflowMode protein_mutation.")
    elif args['alignmentsYAMLinFile'] is None:
        #check that alignment atoms for both ligands are given if one is given
        if args['LIG1refatoms'] is not None and args['LIG2refatoms'] is None:
            parser.error("--LIG1refatoms requires --LIG2refatoms to be specified as well.")
        #convert reference atoms to lists
        lig1_ref_atoms = [int(i) for i in args['LIG1refatoms'].split(',')]
        lig2_ref_atoms = [int(i) for i in args['LIG2refatoms'].split(',')]
    else:
        with open(args['alignmentsYAMLinFile'], 'r')  as f:
            alignments = yaml.safe_load(f)

    run_atm(options, jobname, receptor_file, lig1_file, lig2_file,
            workflow_mode = workflow_mode,
            alignments = alignments,
            pair_alignments = pair_alignments,
            lig1_ref_atoms = lig1_ref_atoms, lig2_ref_atoms = lig2_ref_atoms,
            ff_json_file = ff_json_file)
