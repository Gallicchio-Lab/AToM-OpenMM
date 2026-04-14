import os
import argparse
import yaml
import string
import pickle
import numpy as np
from pathlib import Path
import pandas as pd
from rdkit import Chem
import openmm as mm
from openmm import *
from openmm.app import *
from openmm.unit import *
#from multiprocessing import freeze_support
from atom_openmm.make_atm_system_from_rcpt_lig import make_system
from atom_openmm.rbfe_structprep import rbfe_structprep
from atom_openmm.rbfe_production import rbfe_production
from atom_openmm.utils.AtomUtils import AtomUtils, get_selected_principal_groups
from atom_openmm.uwham import calculate_uwham_from_rundir, create_quality_assessment_plot 
from make_pp_indexes import make_pp_indexes
#from openmm_run import openmm_run

def get_indexes_from_query(topology, query):
    indexes = [ atom.index for atom in topology.atoms() if eval(query, {"atom": atom}) ]
    indexes.sort()
    return indexes

#get the indexes of the atoms of a residue. Optionally, filter them by a query
def get_indexes_from_residue(residue, query = 'True'):
    indexes = [ atom.index for atom in residue.atoms() if eval(query, {"atom": atom}) ]
    indexes.sort()
    return indexes

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

#coordinates of a pdb file without relying on RDKit chemistry perception
def get_pdb_coords(solute_fpath: Path):
    coords = []
    with open(solute_fpath, 'r') as f:
        for line in f:
            if line.startswith(('ATOM  ', 'HETATM')):
                coords.append([
                    float(line[30:38]),
                    float(line[38:46]),
                    float(line[46:54])
                ])
    return np.array(coords)

def _get_solute_coords(solute_fpath: Path):  # pdb or sdf file format
    """
    Return N*3 array for solute coordinates
    """
    assert solute_fpath.suffix in [".sdf", ".pdb"], f"{solute_fpath} is not supported."
    if solute_fpath.suffix == ".sdf":
        mol = Chem.SDMolSupplier(str(solute_fpath), removeHs=False)[0]

        conf = mol.GetConformer()
        N_atoms = mol.GetNumAtoms()
        coords = np.zeros((N_atoms, 3))
        for row in range(N_atoms):
            coords[row] = np.array(list(conf.GetAtomPosition(row)))

        return coords
    else:
        return get_pdb_coords(solute_fpath)

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

    #CM atoms of the receptor, all C-alpha atoms except the protein partners
    rcpt_frame_query = 'atom.residue.chain.id not in ["L","M"] and atom.name == "CA"'
    rcpt_frame_indexes = get_indexes_from_query(topology, rcpt_frame_query)

    #internal receptor frame
    rcpt_frame = get_selected_principal_groups(topology, positions, rcpt_frame_indexes)

    #protein-protein ligand indexes from the generated system topology
    pp_indexes = make_pp_indexes(
        topology = topology,
        chainLig1 = 'L',
        chainLig2 = 'M',
        residLig1 = options['MUTATION_RESID'],
        residLig2 = options['MUTATION_RESID'],
        rcpt_cm_atoms = rcpt_frame['origin']['indices'],
        pos_restrained_atoms = rcpt_frame['origin']['indices'])
    options.update(pp_indexes)

    #LIGAND1_ATOMS (partner chain L)
    ligand1_atoms = options['LIGAND1_ATOMS']

    #LIGAND2_ATOMS (partner chain M)
    ligand2_atoms = options['LIGAND2_ATOMS']

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
        topology, '(atom.residue.chain.id == "B") and (atom.element.atomic_number != 1)')
    options['EXCLUSION_POT_MOL2_INDEXES'] = get_indexes_from_query(
        topology, '(atom.residue.chain.id == "M") and (atom.element.atomic_number != 1)')

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
            alignments = None,
            lig1_ref_atoms = None,
            lig2_ref_atoms = None,
            ff_json_file = 'ff.json',
            csv_datafileout_leg1 = None,
            csv_datafileout_leg2 = None,
            figfileout = None,
            mutation_residue = None):

    #basename
    options['BASENAME'] = jobname

    #protein-protein mutation residue used to define the common/variable split
    options['MUTATION_RESID'] = str(mutation_residue)

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
        ddg, ddg_std, uwham_data = calculate_uwham_from_rundir(
            options['WORKDIR'], options['BASENAME'], mintimeid = discard)

        dgbind1 = uwham_data['dg_leg1']
        dgbind1_std = uwham_data['dg_stderr_leg1']
        dgbind2 = uwham_data['dg_leg2']
        dgbind2_std = uwham_data['dg_stderr_leg2']
        nsamples = uwham_data['nsamples']
        print(f'Relative binding free energy estimate after {nsamples+discard} perturbation energy samples per replica discarding the first {discard} samples:')
        print(f"{jobname}: DG = {ddg:8.3f} +/- {ddg_std:8.3f}  DG(leg1) = {dgbind1:8.3f} +/- {dgbind1_std:8.3f}   DG(leg2) = {dgbind2:8.3f} +/- {dgbind2_std:8.3f} kcal/mol, n_samples: {nsamples}")

        #creates a plot for quality assessment
        df1 = uwham_data['df_leg1']
        N = len(uwham_data['uwham_out_leg1']['W'][:,0])
        df1['W'] = uwham_data['uwham_out_leg1']['W'][:,0]/float(N)
        df2 = uwham_data['df_leg2']
        N = len(uwham_data['uwham_out_leg2']['W'][:,0])
        df2['W'] = uwham_data['uwham_out_leg2']['W'][:,0]/float(N)
        if csv_datafileout_leg1:
            df1.to_csv(csv_datafileout_leg1, index=False)
        if csv_datafileout_leg2:
            df2.to_csv(csv_datafileout_leg2, index=False)
        if figfileout:
            fig = create_quality_assessment_plot(df1, df2)
            fig.savefig(figfileout)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    
    #required arguments
    parser.add_argument('--optionsYAMLinFile',        help='A YAML file with the default options (alchemical schedule, etc.)', required=True)
    parser.add_argument('--jobBasename',              help='The basename of the job', required=True)
    parser.add_argument('--receptorinFile',           help='The structure file of the receptor usually in PDB format', required=True)
    parser.add_argument('--LIG1inFile',               help='The structure file of the first protein partner usually in PDB format', required=True)
    parser.add_argument('--LIG2inFile',               help='The structure file of the second protein partner usually in PDB format', required=True)
    parser.add_argument('--mutationResid',            help='The residue id of the shared mutation site in the protein partners', required=True)

    #optional arguments
    parser.add_argument('--forcefieldJSONCachefile', help='file to store force field information', default='ff.json')
    parser.add_argument('--leg1DataCSVoutFile', type=str,
                        help='Saves the samples for leg1 that were processed in a CSV file, includes the WHAM weights')
    parser.add_argument('--leg2DataCSVoutFile', type=str,
                        help='Saves the samples for leg2 that were processed in a CSV file, includes the WHAM weights')
    parser.add_argument('--plotOutFile', type=str,
                        help='A png file where to save the plot for simulation quality assessment')

    #dictionary of arguments
    args = vars(parser.parse_args())

    jobname = args['jobBasename']
    lig1_file = args['LIG1inFile']
    lig2_file = args['LIG2inFile']
    receptor_file = args['receptorinFile']
    ff_json_file = args['forcefieldJSONCachefile']
    
    #load default settings
    options = None
    with open(args['optionsYAMLinFile'], 'r') as file:
        options = yaml.safe_load(file)
    assert(options)

    run_atm(options, jobname, receptor_file, lig1_file, lig2_file,
            ff_json_file = ff_json_file,
            csv_datafileout_leg1 = args['leg1DataCSVoutFile'],
            csv_datafileout_leg2 = args['leg2DataCSVoutFile'],
            figfileout = args['plotOutFile'],
            mutation_residue = args['mutationResid'])
