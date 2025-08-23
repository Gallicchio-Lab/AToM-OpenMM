#! python

# Usage: python make_atm_rbfe_system_frompdb.py <options>
# Emilio Gallicchio, 5/2023 adapted from code by Bill Swope, 11/2021

############################################
#                                          #
#   IMPORTS                                #
#                                          #
############################################

import os, sys
from datetime import datetime
from time import time

# following for argument passing tools
import argparse

# OpenMM components
from openmm import XmlSerializer
from openmm import Vec3
from openmm.app import PDBFile
from openmm.app import ForceField, Modeller
from openmm.app import PME, HBonds, NoCutoff

# OpenFF components from the toolkit
from openff.toolkit.topology import Molecule

from openmm.unit import angstrom, nanometer, amu, molar

# OpenFF and OpenMM components for ligand force field parameters
from openff.toolkit.topology import Molecule


def boundingBoxSizes(positions):
    xmin = positions[0][0]
    xmax = positions[0][0]
    ymin = positions[0][1]
    ymax = positions[0][1]
    zmin = positions[0][2]
    zmax = positions[0][2]
    for i in range(len(positions)):
        x = positions[i][0]
        y = positions[i][1]
        z = positions[i][2]
        # print('type of x ', type(x))
        # print('Site ', i, ' Coord ', positions[i])
        if(x > xmax):
            xmax = x
        if(x < xmin):
            xmin = x
        if(y > ymax):
            ymax = y
        if(y < ymin):
            ymin = y
        if(z > zmax):
            zmax = z
        if(z < zmin):
            zmin = z
    return [ (xmin,xmax), (ymin,ymax), (zmin,zmax) ] 


def make_system(
        receptorfile,
        displacement,
        xmloutfile,
        pdboutfile,
        lig1file=None,
        lig2file=None,
        lig1sdffile=None,
        lig2sdffile=None,
        cofsdffile=None,
        proteinforcefield='amber14-all.xml',
        solventforcefield='amber14/tip3p.xml',
        ligandforcefield='openff-2.0.0',
        ffcachefile=None,
        implsolv=None,
        hmass=1.0,
        ionicstrength=0.15,
        flagverbose=False
    ):
    print('Generate ATM RBFE OpenMM System')
    today = datetime.today()
    print('\nDate and time at start: ', today.strftime('%c'))
    program_start_timer = time()

    if lig1sdffile is not None:
        print('Warning: LIG1SDFinFile id deprecated. Use LIG1inFile')
        if lig1file is None:
            lig1file = lig1sdffile

    if lig2sdffile is not None:
        print('Warning: LIG2SDFinFile id deprecated. Use LIG2inFile')
        if lig2file is None:
            lig2file = lig2sdffile
            
    #catch abfe or rbfe
    rbfe = False
    if lig2file  is not None:
        rbfe = True

    if isinstance(displacement, str):
        displacement = [float(r) for r in displacement.split()]
    displacement = Vec3(*displacement) * angstrom

    #implicit solvent
    if implsolv == 'None':
        implsolv = None

    hmass = float(hmass)


    #####################################################
    #   Echo out the suppliable input parameters        #
    #####################################################


    print('\nUser-supplied input parameters')
    print('Receptor file name:                 ', receptorfile)
    print('Protein force field:                ', proteinforcefield)
    print('Solvent/ion force field             ', solventforcefield )
    print('Ligand force field:                 ', ligandforcefield)
    print('Ligand 1 file name:                 ', lig1sdffile)
    if rbfe:
        print('Ligand 2 file name:                 ', lig2sdffile)
        print('Displacement                        ', displacement)
    print('Topology PDB output file:           ', pdboutfile)
    print('System XML output file:             ', xmloutfile)
    print('Force field cache file:             ', ffcachefile)


    print('Call ForceField for protein and water')
    forcefield = ForceField(proteinforcefield,solventforcefield)
    if implsolv is not None:
        if implsolv == "OBC2":
            forcefield.loadFile('implicit/obc2.xml')
        elif implsolv == "GBN2":
            forcefield.loadFile('implicit/gbn2.xml')
        elif implsolv == "HCT":
            forcefield.loadFile('implicit/hct.xml')
        elif implsolv == "Vacuum" or implsolv == "vacuum":
            pass
        else:
            print('Unknown implicit solvent %s' % implsolv)
            sys.exit(1)

    # to store OpenFF molecule objects of non-protein units
    ligandmolecules = []

    ############################################
    #                                          #
    #   READ AND CHARACTERIZE RECEPTOR         #
    #                                          #
    ############################################

    rcptpext = os.path.splitext(receptorfile)[1]
    rcpt_ommtopology = None
    rcpt_positions = None
    if rcptpext == '.pdb':
        print('Receptor in PDB format')
        pdbrcpt = PDBFile(receptorfile)
        rcpt_positions = pdbrcpt.positions
        rcpt_ommtopology = pdbrcpt.topology
    elif rcptpext == '.sdf':
        print('Receptor in SDF format')
        molrcpt = Molecule.from_file(receptorfile, file_format='SDF',
                                    allow_undefined_stereo=True)
        ligandmolecules.append(molrcpt)

        pos = molrcpt.conformers[0].to('angstrom').magnitude
        rcpt_positions = [Vec3(pos[i][0], pos[i][1], pos[i][2]) for i in range(pos.shape[0])] * angstrom
        
        offtopology = molrcpt.to_topology()
        rcpt_ommtopology = offtopology.to_openmm(ensure_unique_atom_names=True)
    else:
        print("Error: Unrecognized receptor file name: %s" % receptorfile)
        sys.exit(1)

    nrcpt = rcpt_ommtopology.getNumAtoms()
    print('Number of atoms in receptor:', nrcpt)

    print('Call Modeller: include receptor')
    modeller = Modeller(rcpt_ommtopology, rcpt_positions)

    print("Calculating receptor bounding box:")
    bbox = boundingBoxSizes(rcpt_positions)
    bboxsizes = [ bbox[i][1]-bbox[i][0] for i in range(3) ]
    bboxfaces = [ bboxsizes[2]*bboxsizes[1], bboxsizes[2]*bboxsizes[0],  bboxsizes[1]*bboxsizes[0] ]
    print("Areas of faces", bboxfaces)
    smallest_direction = 0
    smallest_area = bboxfaces[0]
    for i in range(3):
        if bboxfaces[i] < smallest_area:
            smallest_direction = i
    print("Direction of smallest area dimension:", smallest_direction)

        
    ############################################
    #                                          #
    #   READ AND CHARACTERIZE LIGANDS          #
    #                                          #
    ############################################

    if cofsdffile is not None:
        print('Read cofactor:')
        molcof = Molecule.from_file(cofsdffile, file_format='SDF',
                                allow_undefined_stereo=True)
        ligandmolecules.append(molcof)
        molcof_ommtopology = molcof.to_topology().to_openmm(ensure_unique_atom_names=True)

        #assign the residue name, assumes one residue
        resfile = os.path.split(cofsdffile)[1]
        resname = os.path.splitext(resfile)[0]
        for residue in molcof_ommtopology.residues():
            residue.name = resname.upper()

        pos = molcof.conformers[0].to('angstrom').magnitude
        molcof_positions = [Vec3(pos[i][0], pos[i][1], pos[i][2]) for i in range(pos.shape[0])] * angstrom
        ncof = molcof_ommtopology.getNumAtoms()
        print('Number of atoms in cofactor:', ncof)
        print('Call Modeller: include cofactor')
        modeller.add(molcof_ommtopology, molcof_positions)


    print('Read ligand 1:')
    fileext = (os.path.splitext(lig1file)[1]).upper()
    if fileext == '.SDF':
        mollig1 = Molecule.from_file(lig1file, file_format='SDF', allow_undefined_stereo=True)
        ligandmolecules.append(mollig1)
        lig1_ommtopology = mollig1.to_topology().to_openmm(ensure_unique_atom_names=True)
        pos = mollig1.conformers[0].to('angstrom').magnitude
        lig1_positions = [Vec3(pos[i][0], pos[i][1], pos[i][2]) for i in range(pos.shape[0])] * angstrom
        #assign the residue name, assumes one residue
        resname_lig1 = "L1"
        for residue in lig1_ommtopology.residues():
            residue.name = resname_lig1
    elif fileext == '.PDB':
        lig1pdb = PDBFile(lig1file)
        lig1_ommtopology = lig1pdb.topology
        lig1_positions = lig1pdb.positions
        chainname_lig1 = "L"
        for chain in lig1_ommtopology.chains():
            chain.id = chainname_lig1
    else:
        print("Error: unrecognized file: %s" % lig1file)
        sys.exit(1)

    nlig1 = lig1_ommtopology.getNumAtoms()
    print('Number of atoms in ligand 1:', nlig1)
    print('Call Modeller: include ligand 1')
    modeller.add(lig1_ommtopology, lig1_positions)

    if not rbfe:
        # if ABFE, translate the ligand 1 coordinates into the solvent to calculate the
        # bounding box below
        for i in range(nlig1):
            lig1_positions[i] += displacement
    else:
        # RBFE mode:
        # read ligand 2 and place it in the solvent
        print('Read ligand 2:')
        fileext = (os.path.splitext(lig2file)[1]).upper()
        if fileext == '.SDF':
            mollig2 = Molecule.from_file(lig2sdffile, file_format='SDF',
                                         allow_undefined_stereo=True)
            ligandmolecules.append(mollig2)
            lig2_ommtopology = mollig2.to_topology().to_openmm(ensure_unique_atom_names=True)
            pos = mollig2.conformers[0].to('angstrom').magnitude
            lig2_positions = [Vec3(pos[i][0], pos[i][1], pos[i][2]) for i in range(pos.shape[0])] * angstrom
            #assign the residue name, assumes one residue
            resname_lig2 = "L2"
            for residue in lig2_ommtopology.residues():
                residue.name = resname_lig2
        elif fileext == '.PDB':
            lig2pdb = PDBFile(lig2file)
            lig2_ommtopology = lig2pdb.topology
            lig2_positions = lig2pdb.positions
            chainname_lig2 = "M"
            for chain in lig2_ommtopology.chains():
                chain.id = chainname_lig2
        else:
            print("Error: unrecognized file: %s" % lig2file)
            sys.exit(1)
        
        nlig2 = lig2_ommtopology.getNumAtoms()
        print('Number of atoms in ligand 2:', nlig2)
        for i in range(nlig2):
            lig2_positions[i] += displacement
        print('Call Modeller: include ligand 2')
        modeller.add(lig2_ommtopology, lig2_positions)
        lig2atom_indexes = [ i for i in range(nrcpt+nlig1,nrcpt+nlig1+nlig2)]

    print("Calculating system bounding box:")
    if not rbfe:
        bbox = boundingBoxSizes(rcpt_positions + lig1_positions)
    else:
        bbox = boundingBoxSizes(rcpt_positions + lig1_positions + lig2_positions)
    bboxsizes = [ bbox[i][1]-bbox[i][0] for i in range(3) ]
    padding = 2. * 1.0*nanometer
    xBoxvec = Vec3((bboxsizes[0]+padding)/nanometer, 0., 0.)*nanometer
    yBoxvec = Vec3(0.0, (bboxsizes[1]+padding)/nanometer, 0.)*nanometer
    zBoxvec = Vec3(0.0, 0.0, (bboxsizes[2]+padding)/nanometer)*nanometer
    print("boxVectors:", (xBoxvec,yBoxvec,zBoxvec ))



    #bboxfaces = [ bboxsizes[2]*bboxsizes[1], bboxsizes[2]*bboxsizes[0],  bboxsizes[1]*bboxsizes[0] ]
    #print("Areas of faces", bboxfaces)
    #smallest_direction = 0
    #smallest_area = bboxfaces[0]
    #for i in range(3):
    #    if bboxfaces[i] < smallest_area:
    #        smallest_direction = i
    #print("Smallest direction", smallest_direction)


        
    ############################################
    #                                          #
    #   SET UP FORCEFIELD FOR LIGANDS          #
    #                                          #
    ############################################

    print('\nSet up the combined protein + ligand + water system for simulation')
    template_gen = None
    if ligandforcefield[0:4] == "gaff":
        from openmmforcefields.generators import GAFFTemplateGenerator
        print('Using GAFFTemplateGenerator function for ligands')
        template_gen = GAFFTemplateGenerator(molecules=ligandmolecules, cache=ffcachefile )
    elif ligandforcefield[0:6] == "openff":
        from openmmforcefields.generators import SMIRNOFFTemplateGenerator
        print('Call SMIRNOFFTemplateGenerator function for ligands')
        template_gen = SMIRNOFFTemplateGenerator(molecules=ligandmolecules, forcefield=ligandforcefield, cache=ffcachefile )
    elif ligandforcefield[0:8] == "espaloma":
        from openmmforcefields.generators import EspalomaTemplateGenerator
        print('Call EspalomaTemplateGenerator function for ligands')
        template_gen = EspalomaTemplateGenerator(molecules=ligandmolecules, forcefield=ligandforcefield, cache=ffcachefile )
    else:
        print('Unknown ligand force field %s' % ligandforcefield)
        sys.exit(1)

    # Register the SMIRNOFF template generator
    # NOTE: forcefield object was initialized (above)
    # for the protein + water (subsystem) force field 
    # This step adds support for the ligand force field
    forcefield.registerTemplateGenerator(template_gen.generator)

    if implsolv is None:
        print("Ionic strength = ", ionicstrength*molar)
        print("Adding solvent and processing system ...")
        modeller.addSolvent(forcefield, boxVectors = (xBoxvec,yBoxvec,zBoxvec ), ionicStrength = ionicstrength*molar)
        print("Number of atoms in solvated system:", modeller.topology.getNumAtoms())
        system=forcefield.createSystem(modeller.topology, nonbondedMethod = PME, nonbondedCutoff = 0.9*nanometer,
                                    constraints=HBonds, rigidWater = True, removeCMMotion = False, hydrogenMass = hmass*amu)
    else:
        print("Solvent model: %s" % implsolv)
        print("Number of atoms in implicit solvent system:", modeller.topology.getNumAtoms())
        print("Processing system ...")
        system=forcefield.createSystem(modeller.topology, nonbondedMethod = NoCutoff,
                                    constraints=HBonds, rigidWater = True, removeCMMotion = False, hydrogenMass = hmass*amu)

    with open(xmloutfile, 'w') as output:
        output.write(XmlSerializer.serialize(system))

    if pdboutfile is not None:
        PDBFile.writeFile(modeller.topology, modeller.positions,
                        open(pdboutfile,'w'), keepIds=True)

    today = datetime.today()
    print('\n\nDate and time at end:   ', today)
    program_end_timer = time()
    print('\nTotal compute time %.3f seconds' %  (program_end_timer-program_start_timer))


def main():
    whatItDoes = """
    Produces an .xml file with the OpenMM's system for ATM relative binding
    free energy calcuations.

    The user supplies a pdb file that contains protein, cofactors, ligands, 
    water, ions, probably prepared by a different software package 
    (Maestro, OpenEye).  Cofactors and ligands have to be described in an
    sdf file.  Then force fields are assigned to all components

    Emilio Gallicchio 5/2023
    adapted from simProteinLigandWater.py script by Bill Swope 11/2021
    """

    parser = argparse.ArgumentParser(description=whatItDoes)

    # Required input
    parser.add_argument('--receptorinFile', required=True,  type=str, default=None, dest='receptorfile',
                        help='Receptor in SDF (.sdf) or PDB format (.pdb)')
    parser.add_argument('--displacement',  required=True,  type=str, default=None, dest='displacement',
                        help='string with displacement vector in angstroms like "22. 0.0 0.0" ')
    parser.add_argument('--systemXMLoutFile',  required=True,  type=str, default=None, dest='xmloutfile',
                        help='Name of the XML file where to save the System')
    parser.add_argument('--systemPDBoutFile', required=True, type=str, default=None, dest='pdboutfile',
                        help='Name of the PDB file where to output to system')

    # Optional input
    parser.add_argument('--LIG1inFile',  required=False,  type=str, default=None, dest='lig1file',
                        help='First ligand in SDF (.sdf) or PDB format (.pdb)')
    parser.add_argument('--LIG2inFile',  required=False,  type=str, default=None, dest='lig2file',
                        help='Second ligand in SDF (.sdf) or PDB format (.pdb)')
    parser.add_argument('--LIG1SDFinFile',  required=False,  type=str, default=None, dest='lig1sdffile',
                        help='SDF file of first ligand')
    parser.add_argument('--LIG2SDFinFile',  required=False,  type=str, default=None, dest='lig2sdffile',
                        help='SDF file of second ligand')
    parser.add_argument('--cofactorsSDFFile',  required=False,  type=str, default=None, dest='cofsdffile',
                        help='SDF file with receptor cofactors')

    parser.add_argument('--proteinForceField', required=False, type=str, dest='proteinforcefield',
                        default='amber14-all.xml',
                        help='Force field for protein` amber14-all.xml ')
    parser.add_argument('--solventForceField', required=False, type=str, dest='solventforcefield',
                        default='amber14/tip3p.xml',
                        help='Force field for solvent/ions` amber14/tip3p.xml ')
    parser.add_argument('--ligandForceField', required=False, type=str, dest='ligandforcefield',
                        default='openff-2.0.0',
                        help='Force field for ligand:  openff-2.0.0, gaff, or espaloma-0.3.2')
    parser.add_argument('--implicitSolvent', required=False, type=str, dest='implsolv',
                        default=None,
                        help='Implicit solvent to use: HCT OBC2 GBn2, vacuum.')
    parser.add_argument('--forcefieldJSONCachefile', required=False, type=str, dest='ffcachefile',
                        default=None,
                        help='Force field ligand cache database')
    parser.add_argument('--hmass', required=False, type=float, dest='hmass',
                        default=1.0,
                        help='Hydrogen mass, set it to 1.5 amu to use a 4 fs time-step')
    parser.add_argument('--ionicStrength', required=False, type=float, dest='ionicstrength',
                        default=0.15,
                        help='Total concentration of monoatomic ions to add')
    
    # Arguments that are flags
    parser.add_argument('--verbose', required=False, action='store_true', dest='flagverbose',
                        help='Get more output with this flag')

    args = vars(parser.parse_args())
    make_system(**args)


if __name__ == "__main__":
    main()
