from __future__ import print_function
from __future__ import division
import sys
import time
import math
import random
import logging
import signal
import shutil
import random

from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from datetime import datetime

from openmm_async_re import *
from ommsystem import *

class OMMSystemAmberABFE_zrestr(OMMSystemAmberABFE):
    
    def set_vsite_restraints(self):
        #prmtop is an instance variable with "self.prmtop"; so "prmtop.topology.atoms"
        #in input string from control file need to be replaced with "self.prmtop.topology.atoms"
        #to set the indexes for the receptor and ligand CM atoms
        
        #add receptor and ligand atom indexes for centroid calculation using string based python syntaxes from the command file
        
        cm_lig_atoms_selection = self.keywords.get('LIGAND_CM_ATOMS')   #python syntax as string for selecting specific receptor atoms
        if cm_lig_atoms_selection is not None:
            cm_lig_atoms_selection = cm_lig_atoms_selection.replace("prmtop.topology.atoms", "self.prmtop.topology.atoms")
            #cm_lig_atom_str ="lig_atom_restr = [{0}]".format(cm_lig_atoms_selection)
            cm_lig_atom_str ="[{0}]".format(cm_lig_atoms_selection)
            #print(cm_lig_atom_str)
            lig_atom_restr = eval(cm_lig_atom_str)
            print(lig_atom_restr)
        else:
            lig_atom_restr = None

        cm_rcpt_atoms_selection = self.keywords.get('RCPT_CM_ATOMS')   #python syntax as string for selecting specific ligand atoms
        if cm_rcpt_atoms_selection is not None:
            cm_rcpt_atoms_selection = cm_rcpt_atoms_selection.replace("prmtop.topology.atoms", "self.prmtop.topology.atoms")
            cm_rcpt_atom_str = "[{0}]".format(cm_rcpt_atoms_selection)
            #print(cm_rcpt_atom_str)
            rcpt_atom_restr = eval(cm_rcpt_atom_str)
            print(rcpt_atom_restr)
        else:
            rcpt_atom_restr = None

        cmrestraints_present = (cm_rcpt_atoms_selection is not None) and (cm_lig_atoms_selection is not None)

        self.vsiterestraintForce = None
        if cmrestraints_present:
            print("Adding z-restraints")
            cmkf = float(self.keywords.get('CM_KF'))
            kf = cmkf * kilocalories_per_mole/angstrom**2 #force constant for Vsite CM-CM restraint
            cmtol = float(self.keywords.get('CM_TOL'))
            r0 = cmtol * angstrom #radius of Vsite sphere
            ligoffset = self.keywords.get('LIGOFFSET')
            offz = None
            if ligoffset:
                ligoffset = [float(offset) for offset in ligoffset.split(',')]*angstrom
                offz = ligoffset[2]
            #functional form of tre Z-restraint force
            self.vsiterestraintForce = mm.CustomCentroidBondForce(2,"0.5*kf*( step(d)*max(0,d-r0)^2 + step(-d)*max(0,-d-r0)^2 ) ; d = z2 - offz - z1")
            self.vsiterestraintForce.addPerBondParameter("kf")
            self.vsiterestraintForce.addPerBondParameter("r0")
            self.vsiterestraintForce.addPerBondParameter("offz")
            self.vsiterestraintForce.setForceGroup(1)
            self.vsiterestraintForce.addGroup(rcpt_atom_restr)
            self.vsiterestraintForce.addGroup(lig_atom_restr)
            self.vsiterestraintForce.addBond([0,1], [kf, r0, offz ])
            self.system.addForce(self.vsiterestraintForce)

    def set_ligand_atoms(self):
        #prmtop is an instance variable with "self.prmtop"; so "prmtop.topology.atoms"
        #in input string from control file need to be replaced with "self.prmtop.topology.atoms"
        #to set the indexes for the receptor and ligand CM atoms
        
        lig_atoms_selection = self.keywords.get('LIGAND_ATOMS')
        if lig_atoms_selection is not None:
            lig_atoms_selection = lig_atoms_selection.replace("prmtop.topology.atoms", "self.prmtop.topology.atoms")
            lig_atom_str = "[{0}]".format(lig_atoms_selection)
            self.lig_atoms = eval(lig_atom_str)
            print(self.lig_atoms)
        else:
            msg = "Error: LIGAND_ATOMS is required"
            self._exit(msg)
    
    
class openmm_job_AmberABFE_zrestr(openmm_job_ATM):
    def __init__(self, command_file, options):
        super().__init__(command_file, options)

        prmtopfile = self.basename + ".prmtop"
        crdfile = self.basename + ".inpcrd"

        if not self.stateparams:
            self._buildStates()
        
        #builds service worker for replicas use
        service_ommsys = OMMSystemAmberABFE_zrestr(self.basename, self.keywords, prmtopfile, crdfile, self.logger)
        self.service_worker = OMMWorkerATM(self.basename, service_ommsys, self.keywords, compute = False)
        #creates openmm replica objects
        self.openmm_replicas = []
        for i in range(self.nreplicas):
            replica = OMMReplicaATM(i, self.basename, self.service_worker, self.logger)
            if replica.stateid == None:
                replica.set_state(i, self.stateparams[i])#initial setting
            self.openmm_replicas.append(replica)

        # creates openmm workers objects
        self.openmm_workers = []
        pattern = re.compile('(\d+):(\d+)')
        for node in self.compute_nodes:
            slot_id = node["slot_number"]
            matches = pattern.search(slot_id)
            platform_id = int(matches.group(1))
            device_id = int(matches.group(2))
            ommsys = OMMSystemAmberABFE_zrestr(self.basename, self.keywords, prmtopfile, crdfile, self.logger)
            self.openmm_workers.append(OMMWorkerATM(self.basename, ommsys, self.keywords, self.gpu_platform_name, platform_id, device_id))


if __name__ == '__main__':

    # Parse arguments:
    usage = "%prog <ConfigFile>"

    if len(sys.argv) != 2:
        print("Please specify ONE input file")
        sys.exit(1)

    commandFile = sys.argv[1]

    print("")
    print("====================================")
    print("BEDAM Asynchronous Replica Exchange ")
    print("====================================")
    print("")
    print("Started at: " + str(time.asctime()))
    print("Input file:", commandFile)
    print("")
    sys.stdout.flush()

    rx = openmm_job_AmberABFE_zrestr(commandFile, options=None)

    rx.setupJob()

    rx.scheduleJobs()
