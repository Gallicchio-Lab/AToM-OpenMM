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

from openmm_async_re import openmm_job_AmberRBFE
from openmm_async_re import *
from ommsystem import *

class OMMSystemAmberRBFE_zrestr(OMMSystemAmberRBFE):

    def set_vsite_restraints(self):
        #prmtop is an instance variable with "self.prmtop"; so "prmtop.topology.atoms"
        #in input string from control file need to be replaced with "self.prmtop.topology.atoms"
        #to set the indexes for the receptor and ligand CM atoms
        
        #add receptor and ligand atom indexes for centroid calculation using string based python syntaxes from the command file
        
        cm_lig_atoms_selection = self.keywords.get('LIGAND1_CM_ATOMS')   #python syntax as string for selecting specific receptor atoms
        if cm_lig_atoms_selection is not None:
            cm_lig_atoms_selection = cm_lig_atoms_selection.replace("prmtop.topology.atoms", "self.prmtop.topology.atoms")
            #cm_lig_atom_str ="lig_atom_restr = [{0}]".format(cm_lig_atoms_selection)
            cm_lig_atom_str ="[{0}]".format(cm_lig_atoms_selection)
            #print(cm_lig_atom_str)
            lig1_atom_restr = eval(cm_lig_atom_str)
            print(lig1_atom_restr)
        else:
            lig1_atom_restr = None

        cm_lig_atoms_selection = self.keywords.get('LIGAND2_CM_ATOMS')   #python syntax as string for selecting specific receptor atoms
        if cm_lig_atoms_selection is not None:
            cm_lig_atoms_selection = cm_lig_atoms_selection.replace("prmtop.topology.atoms", "self.prmtop.topology.atoms")
            #cm_lig_atom_str ="lig_atom_restr = [{0}]".format(cm_lig_atoms_selection)
            cm_lig_atom_str ="[{0}]".format(cm_lig_atoms_selection)
            #print(cm_lig_atom_str)
            lig2_atom_restr = eval(cm_lig_atom_str)
            print(lig2_atom_restr)
        else:
            lig2_atom_restr = None


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
            
            self.vsiterestraintForce = mm.CustomCentroidBondForce(2,"0.5*kf*( step(d)*max(0,d-r0)^2 + step(-d)*max(0,-d-r0)^2 ) ; d = z2 - offz - z1")
            self.vsiterestraintForce.addPerBondParameter("kf")
            self.vsiterestraintForce.addPerBondParameter("r0")
            self.vsiterestraintForce.addPerBondParameter("offz")
            self.vsiterestraintForce.addGroup(rcpt_atom_restr)
            self.vsiterestraintForce.addGroup(lig1_atom_restr)
            self.vsiterestraintForce.addGroup(lig2_atom_restr)
            self.vsiterestraintForce.addBond([0,1], [kf/(kilocalories_per_mole/angstrom**2), r0/angstrom, self.lig1offset[2]/angstrom])
            self.vsiterestraintForce.addBond([0,2], [kf/(kilocalories_per_mole/angstrom**2), r0/angstrom, self.lig2offset[2]/angstrom])
            self.system.addForce(self.vsiterestraintForce)

    def set_ligand_atoms(self):
        #prmtop is an instance variable with "self.prmtop"; so "prmtop.topology.atoms"
        #in input string from control file need to be replaced with "self.prmtop.topology.atoms"
        #to set the indexes for the receptor and ligand CM atoms
        
        lig1_atoms_in = self.keywords.get('LIGAND1_ATOMS')   #indexes of ligand1 atoms
        lig2_atoms_in = self.keywords.get('LIGAND2_ATOMS')   #indexes of ligand2 atoms
        if lig1_atoms_in is not None:
            lig1_atoms_in = lig1_atoms_in.replace("prmtop.topology.atoms", "self.prmtop.topology.atoms")
            lig1_atom_str = "[{0}]".format(lig1_atoms_in)
            self.lig1_atoms = eval(lig1_atom_str)
            print(self.lig1_atoms)
        else:
            msg = "Error: LIGAND1_ATOMS is required"
            self._exit(msg)
            
        if lig2_atoms_in is not None:
            lig2_atoms_in = lig2_atoms_in.replace("prmtop.topology.atoms", "self.prmtop.topology.atoms")
            lig2_atom_str = "[{0}]".format(lig2_atoms_in)
            self.lig2_atoms = eval(lig2_atom_str)
            print(self.lig2_atoms)
        else:
            msg = "Error: LIGAND2_ATOMS is required"
            self._exit(msg)

class openmm_job_AmberRBFE_zrestr(openmm_job_ATM):
    def __init__(self, command_file, options):
        super().__init__(command_file, options)

        prmtopfile = self.basename + ".prmtop"
        crdfile = self.basename + ".inpcrd"

        if self.stateparams is None:
            self._buildStates()
        
        #builds service worker for replicas use
        service_ommsys = OMMSystemAmberRBFE_zrestr(self.basename, self.keywords, prmtopfile, crdfile, self.logger)
        self.service_worker = OMMWorkerATM(self.basename, service_ommsys, self.keywords, compute = False, logger = self.logger)
        #creates openmm replica objects
        self.openmm_replicas = []
        for i in range(self.nreplicas):
            replica = OMMReplicaATM(i, self.basename, self.service_worker, self.logger)
            if replica.stateid == None:
                replica.set_state(i, self.stateparams[i])#initial setting
            self.openmm_replicas.append(replica)

        # creates openmm context objects
        self.openmm_workers = []
        for node in self.compute_nodes:
            ommsys = OMMSystemAmberRBFE_zrestr(self.basename, self.keywords, prmtopfile, crdfile, self.logger) 
            self.openmm_workers.append(OMMWorkerATM(self.basename, ommsys, self.keywords, node_info = node, compute = True, logger = self.logger))

        
if __name__ == '__main__':

    # Parse arguments:
    usage = "%prog <ConfigFile>"

    if len(sys.argv) != 2:
        print("Please specify ONE input file")
        sys.exit(1)

    commandFile = sys.argv[1]

    print("")
    print("=================================================")
    print("AToM Membrane RTFE Asynchronous Replica Exchange ")
    print("=================================================")
    print("")
    print("Started at: " + str(time.asctime()))
    print("Input file:", commandFile)
    print("")
    sys.stdout.flush()

    rx = openmm_job_AmberRBFE_zrestr(commandFile, options=None)

    rx.setupJob()

    rx.scheduleJobs()
