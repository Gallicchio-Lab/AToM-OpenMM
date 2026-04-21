import os
import sys
import argparse
import yaml
from pathlib import Path
import pandas as pd
import openmm as mm
from openmm import *
from openmm.app import *
from openmm.unit import *
from atom_openmm.make_atm_system_from_rcpt_lig import make_system
from atom_openmm.rbfe_structprep import rbfe_structprep
from atom_openmm.rbfe_production import rbfe_production
from atom_openmm.utils.AtomUtils import (
    calc_displ_vec,
    get_indexes_from_query,
    get_indexes_from_residue,
    get_selected_principal_groups,
)
from atom_openmm.uwham import calculate_uwham_from_rundir, create_quality_assessment_plot

def get_first_residue_by_name(topology, residue_name):
    for chain in topology.chains():
        for residue in chain.residues():
            if residue.name == residue_name:
                return residue
    return None


def get_ligand_geometric_center(positions, residue):
    atom_indexes = get_indexes_from_residue(residue)
    center = Vec3(0.0, 0.0, 0.0) * nanometer
    for atom_index in atom_indexes:
        center += positions[atom_index]
    return center / float(len(atom_indexes))


def squared_distance(position1, position2):
    displacement = (position1 - position2).value_in_unit(nanometer)
    return displacement.x ** 2 + displacement.y ** 2 + displacement.z ** 2


def get_attach_atom_from_residue(residue, positions, attach_index=None):
    atoms = list(residue.atoms())
    if attach_index is not None:
        return atoms[attach_index]
    ligand_center = get_ligand_geometric_center(positions, residue)
    return min(atoms, key=lambda atom: squared_distance(positions[atom.index], ligand_center))


def system_has_ghost_pair(pdb_file):
    pdb = PDBFile(pdb_file)
    residue_names = [residue.name for residue in pdb.topology.residues()]
    return residue_names.count("L1") == 1 and residue_names.count("L2") == 1


def patch_system_with_ghost(pdb_file, xml_file, displacement, ghost_mass, attach_index=None):
    pdb = PDBFile(pdb_file)
    positions_nm = [Vec3(pos.x, pos.y, pos.z) for pos in pdb.positions.value_in_unit(nanometer)]
    topology = pdb.topology

    with open(xml_file) as input:
        system = XmlSerializer.deserialize(input.read())

    ligand_residue = get_first_residue_by_name(topology, "L1")
    assert ligand_residue is not None, "Could not find ligand residue L1 in the prepared system."
    ligand_attach_atom = get_attach_atom_from_residue(
        ligand_residue,
        pdb.positions,
        attach_index=attach_index,
    )
    ligand_atom_indexes = get_indexes_from_residue(ligand_residue)

    bound_anchor_position = Vec3(
        positions_nm[ligand_attach_atom.index].x,
        positions_nm[ligand_attach_atom.index].y,
        positions_nm[ligand_attach_atom.index].z,
    )
    displacement_vec = (Vec3(*displacement) * angstrom).value_in_unit(nanometer)

    for atom_index in ligand_atom_indexes:
        positions_nm[atom_index] += displacement_vec

    ligand_residue.name = "L2"
    ligand_residue.chain.id = "M"

    ghost_topology = Topology()
    ghost_chain = ghost_topology.addChain(id="L")
    ghost_residue = ghost_topology.addResidue("L1", ghost_chain)
    ghost_element = Element.getBySymbol("C")
    ghost_topology.addAtom("C1", ghost_element, ghost_residue)

    modeller = Modeller(topology, positions_nm * nanometer)
    modeller.add(ghost_topology, [bound_anchor_position] * nanometer)

    system.addParticle(ghost_mass * amu)
    for force in system.getForces():
        if isinstance(force, mm.NonbondedForce):
            force.addParticle(0.0, 1.0 * angstrom, 0.0 * kilojoule_per_mole)
        elif isinstance(force, mm.CustomNonbondedForce):
            force.addParticle([0.0] * force.getNumPerParticleParameters())
        elif isinstance(force, mm.GBSAOBCForce):
            force.addParticle(0.0, 1.0, 0.0)

    with open(xml_file, "w") as output:
        output.write(XmlSerializer.serialize(system))
    with open(pdb_file, "w") as output:
        PDBFile.writeFile(modeller.topology, modeller.positions, output, keepIds=True)


def rbfe_setup(receptor_file, lig_file, ff_json_file, options):
    basename = options['BASENAME']
    setup = {}
    setup['receptorfile'] = receptor_file
    setup['lig1file'] = lig_file
    setup['displacement'] = options['DISPLACEMENT']
    setup['xmloutfile'] = basename + "_sys.xml"
    setup['pdboutfile'] = basename + ".pdb"
    setup['ffcachefile'] = ff_json_file
    if 'HMASS' in options:
        setup['hmass'] = options['HMASS']
    if 'LIGAND_FORCE_FIELD' in options:
        setup['ligandforcefield'] = options['LIGAND_FORCE_FIELD']

    make_system(**setup)
    patch_system_with_ghost(
        setup['pdboutfile'],
        setup['xmloutfile'],
        options['DISPLACEMENT'],
        options.get('GHOST_MASS', 12.011),
        attach_index=options.get('LIGAND_ATTACH_INDEX'),
    )


def rbfe_prepare_args(options):
    basename = options['BASENAME']

    pdb = PDBFile(basename + ".pdb")
    topology = pdb.topology
    positions = pdb.positions

    ghost_residue = get_first_residue_by_name(topology, "L1")
    ligand_residue = get_first_residue_by_name(topology, "L2")
    assert ghost_residue is not None
    assert ligand_residue is not None

    ligand1_atoms = get_indexes_from_residue(ghost_residue)
    ligand2_atoms = get_indexes_from_residue(ligand_residue)
    options['LIGAND1_ATOMS'] = ligand1_atoms
    options['LIGAND2_ATOMS'] = ligand2_atoms

    options['LIGAND1_VAR_ATOMS'] = options['LIGAND1_ATOMS']
    options['LIGAND2_VAR_ATOMS'] = options['LIGAND2_ATOMS']

    ghost_attach_atom = list(ghost_residue.atoms())[0]
    ligand_attach_atom = get_attach_atom_from_residue(
        ligand_residue,
        positions,
        attach_index=options.get('LIGAND_ATTACH_INDEX'),
    )
    options['LIGAND1_ATTACH_ATOM'] = ghost_attach_atom.index
    options['LIGAND2_ATTACH_ATOM'] = ligand_attach_atom.index

    options['LIGAND1_CM_ATOMS'] = [ghost_attach_atom.index]
    options['LIGAND2_CM_ATOMS'] = [ligand_attach_atom.index]

    lig1cm_pos = positions[ghost_attach_atom.index]
    lig2cm_pos = positions[ligand_attach_atom.index]
    displ = (lig2cm_pos - lig1cm_pos).value_in_unit(angstrom)
    options['DISPLACEMENT'] = [ displ.x , displ.y, displ.z ] 

    rcpt_chain_name = options.get('RCPT_CHAIN_NAME', 'A')
    rcpt_frame_query = f'atom.residue.chain.id == "{rcpt_chain_name}" and atom.name == "CA"'
    rcpt_frame_indexes = get_indexes_from_query(topology, rcpt_frame_query)
    rcpt_frame = get_selected_principal_groups(topology, positions, rcpt_frame_indexes)

    options['RCPT_CM_ATOMS'] = rcpt_frame['origin']['indices']
    options['RCPT_FRAME_ATOMS_O'] = options['RCPT_CM_ATOMS']
    options['RCPT_FRAME_ATOMS_Z'] = rcpt_frame['z_axis']['indices']
    options['RCPT_FRAME_ATOMS_Y'] = rcpt_frame['y_axis']['indices']

    rcpt_cm_pos = Vec3(
        float(rcpt_frame['origin']['com'][0]),
        float(rcpt_frame['origin']['com'][1]),
        float(rcpt_frame['origin']['com'][2]),
    ) * nanometer
    offset = (lig1cm_pos - rcpt_cm_pos).value_in_unit(angstrom)
    options['LIGOFFSET'] = [ offset.x, offset.y, offset.z ]

    options['POS_RESTRAINED_ATOMS'] = None
    options['EXCLUSION_POT_MOL1_INDEXES'] = get_indexes_from_query(
        topology,
        f'(atom.residue.chain.id == "{rcpt_chain_name}") and (atom.element.atomic_number != 1)',
    )
    options['EXCLUSION_POT_MOL2_INDEXES'] = get_indexes_from_residue(
        ligand_residue, query="(atom.element.atomic_number != 1)"
    )
    options['WORKDIR'] = os.environ.get("PWD")


def rbfe_production_data(options):
    rundir = options['WORKDIR']
    jobname = options['BASENAME']
    nstates = len(options['LAMBDAS'])
    datafiles = [
        os.path.join(rundir, f"r{i}", f"{jobname}.out") for i in range(nstates)
    ]
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
    if line_counts:
        min_file = min(line_counts, key=line_counts.get)
        return line_counts[min_file]
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
    for key, value in replacements.items():
        content = content.replace(key, str(value))
    with open(output_path, 'w') as f:
        f.write(content)

def run_atm(options,
            jobname,
            receptor_file,
            lig_file,
            ff_json_file = 'ff.json',
            csv_datafileout_leg1 = None,
            csv_datafileout_leg2 = None,
            figfileout = None):
    options['BASENAME'] = jobname

    if not options.get('DISPLACEMENT'):
        displ = calc_displ_vec(receptor_file, lig_file)
        options['DISPLACEMENT'] = list(displ)
    assert(options['DISPLACEMENT'])

    needs_setup = (
        not Path(jobname + ".pdb").exists()
        or not Path(jobname + "_sys.xml").exists()
        or not system_has_ghost_pair(jobname + ".pdb")
    )
    if needs_setup:
        rbfe_setup(receptor_file, lig_file, ff_json_file, options)

    rbfe_prepare_args(options)

    if not Path(options['BASENAME'] + ".yaml").exists():
        with open(options['BASENAME'] + ".yaml", "w") as file:
            yaml.dump(options, file, default_flow_style=None, width=1000000, sort_keys=False)

    if Path('vmd_template.in').exists() and not Path('vmd.in').exists():
        create_vmd_infile(options)

    if not Path(options['BASENAME'] + "_0.xml").exists():
        rbfe_structprep(config_file=None, options=options)

    nsamples_per_replica = rbfe_production_data(options)
    if (
        not options.get('MAX_SAMPLES', None)
        or not nsamples_per_replica
        or nsamples_per_replica < options['MAX_SAMPLES']
    ):
        rbfe_production(config_file=None, options=options)

    nsamples_per_replica = rbfe_production_data(options)
    if nsamples_per_replica:
        discard = int(nsamples_per_replica / 3)
        dgb, dgb_std, uwham_data = calculate_uwham_from_rundir(
            options['WORKDIR'], options['BASENAME'], mintimeid=discard
        )

        dg_ghost = uwham_data['dg_leg1']
        dg_ghost_std = uwham_data['dg_stderr_leg1']
        dg_ligand = uwham_data['dg_leg2']
        dg_ligand_std = uwham_data['dg_stderr_leg2']
        nsamples = uwham_data['nsamples']
        print(
            f"Binding free energy estimate after {nsamples + discard} perturbation energy samples "
            f"per replica discarding the first {discard} samples:"
        )
        print(
            f"{jobname}: DGb = {dgb:8.3f} +/- {dgb_std:8.3f}  "
            f"DG(ghost leg) = {dg_ghost:8.3f} +/- {dg_ghost_std:8.3f}  "
            f"DG(ligand leg) = {dg_ligand:8.3f} +/- {dg_ligand_std:8.3f} "
            f"kcal/mol, n_samples: {nsamples}"
        )

        df1 = uwham_data['df_leg1']
        n1 = len(uwham_data['uwham_out_leg1']['W'][:, 0])
        df1['W'] = uwham_data['uwham_out_leg1']['W'][:, 0] / float(n1)
        df2 = uwham_data['df_leg2']
        n2 = len(uwham_data['uwham_out_leg2']['W'][:, 0])
        df2['W'] = uwham_data['uwham_out_leg2']['W'][:, 0] / float(n2)
        if csv_datafileout_leg1:
            df1.to_csv(csv_datafileout_leg1, index=False)
        if csv_datafileout_leg2:
            df2.to_csv(csv_datafileout_leg2, index=False)
        if figfileout:
            fig = create_quality_assessment_plot(df1, df2)
            fig.savefig(figfileout)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--optionsYAMLinFile", required=True)
    parser.add_argument("--jobBasename", required=True)
    parser.add_argument("--receptorinFile", required=True)
    parser.add_argument("--LIG1inFile", required=True)

    parser.add_argument("--forcefieldJSONCachefile", default="ff.json")
    parser.add_argument("--leg1DataCSVoutFile", type=str)
    parser.add_argument("--leg2DataCSVoutFile", type=str)
    parser.add_argument("--plotOutFile", type=str)
    parser.add_argument(
        "--ligandAttachAtomIndex",
        type=int,
        help="Optional 1-based ligand atom index to use as the attachment atom.",
    )

    args = vars(parser.parse_args())

    with open(args['optionsYAMLinFile'], "r") as file:
        options = yaml.safe_load(file)

    if args['ligandAttachAtomIndex'] is not None:
        options['LIGAND_ATTACH_INDEX'] = int(args['ligandAttachAtomIndex']) - 1

    run_atm(
        options,
        args['jobBasename'],
        args['receptorinFile'],
        args['LIG1inFile'],
        ff_json_file=args['forcefieldJSONCachefile'],
        csv_datafileout_leg1=args['leg1DataCSVoutFile'],
        csv_datafileout_leg2=args['leg2DataCSVoutFile'],
        figfileout=args['plotOutFile'],
    )
