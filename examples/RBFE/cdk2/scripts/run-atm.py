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
from atom_openmm.utils.AtomUtils import AtomUtils, get_selected_principal_groups, get_indexes_from_query, get_indexes_from_residue, cm_from_indexes, calc_displ_vec
from atom_openmm.uwham import calculate_uwham_from_rundir, create_quality_assessment_plot
#from openmm_run import openmm_run

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

    #LIGAND1_ATOMS (residue L1)
    res1 = None
    lig1resname = 'L1'
    for chain in topology.chains():
        for residue in chain.residues():
            if residue.name == lig1resname:
                res1 = residue
                break
    assert(res1)
    ligand1_atoms = get_indexes_from_residue(res1)
    options['LIGAND1_ATOMS'] = ligand1_atoms

    #LIGAND2_ATOMS (residue L2)
    res2 = None
    lig2resname = 'L2'
    for chain in topology.chains():
        for residue in chain.residues():
            if residue.name == lig2resname:
                res2 = residue
                break
    assert(res2)
    ligand2_atoms = get_indexes_from_residue(res2)
    options['LIGAND2_ATOMS'] = ligand2_atoms

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
    options['EXCLUSION_POT_MOL2_INDEXES'] = get_indexes_from_residue(res2, query = '(atom.element.atomic_number != 1)')

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
            figfileout = None):

    #basename
    options['BASENAME'] = jobname

    #alignment atoms
    if alignments is not None:
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
        displ = calc_displ_vec(receptor_file, lig2_file)
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
    parser.add_argument('--LIG1inFile',               help='The structure file of the first ligand usually in SDF format', required=True)
    parser.add_argument('--LIG2inFile',               help='The structure file of the first ligand usually in SDF format', required=True)

    #optional arguments.
    #the alignment atoms can be specified explicitly or a pickle file with the alignment dictionary must be given
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument( '--alignmentsYAMLinFile',     help='The yaml file with the alignment dictionary', type=str)
    group.add_argument( '--LIG1refatoms',             help='The alignment atoms of the first ligand, --LIG1refatoms "14,12,11" for example, starting with 1', type=str)
    parser.add_argument('--LIG2refatoms',             help='The alignment atoms of the second ligand, --LIG2refatoms "14,12,12" for example, starting with 1', type=str)

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

    alignments = None
    lig1_ref_atoms = None
    lig2_ref_atoms = None
    if args['alignmentsYAMLinFile'] is None:
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
            alignments = alignments,
            lig1_ref_atoms = lig1_ref_atoms, lig2_ref_atoms = lig2_ref_atoms,
            ff_json_file = ff_json_file,
            csv_datafileout_leg1 = args['leg1DataCSVoutFile'],
            csv_datafileout_leg2 = args['leg2DataCSVoutFile'],
            figfileout = args['plotOutFile'])
