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
    def create_system(self):

        self.prmtop = AmberPrmtopFile(self.prmtopfile)
        self.inpcrd = AmberInpcrdFile(self.crdfile)
        self.system = self.prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                                          constraints=HBonds)
        self.topology = self.prmtop.topology
        self.positions = self.inpcrd.positions
        self.boxvectors = self.inpcrd.boxVectors
        
        atm_utils = ATMMetaForceUtils(self.system)

        lig_atoms = self.keywords.get('LIGAND_ATOMS')   #indexes of ligand atoms
        if lig_atoms:
            lig_atoms = [int(i) for i in lig_atoms]
        else:
            msg = "Error: LIGAND_ATOMS is required"
            self._exit(msg)
        
        cm_lig_atoms = self.keywords.get('REST_LIGAND_CMLIG_ATOMS')   #indexes of ligand atoms for CM-CM Vsite restraint
        if cm_lig_atoms:
            lig_atom_restr = [int(i) for i in cm_lig_atoms]
        else:
            lig_atom_restr = None

        cm_rcpt_atoms = self.keywords.get('REST_LIGAND_CMREC_ATOMS')   #indexes of rcpt atoms for CM-CM Vsite restraint
        if cm_rcpt_atoms:
            rcpt_atom_restr = [int(i) for i in cm_rcpt_atoms]
        else:
            rcpt_atom_restr = None

        cmrestraints_present = (cm_rcpt_atoms is not None) and (cm_lig_atoms is not None)
        
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
            zrestrforce = mm.CustomCentroidBondForce(2,"0.5*kf*( step(d)*max(0,d-r0)^2 + step(-d)*max(0,-d-r0)^2 ) ; d = z2 - offz - z1")
            zrestrforce.addPerBondParameter("kf")
            zrestrforce.addPerBondParameter("r0")
            zrestrforce.addPerBondParameter("offz")
            zrestrforce.setForceGroup(1)
            zrestrforce.addGroup(rcpt_atom_restr)
            zrestrforce.addGroup(lig_atom_restr)
            zrestrforce.addBond([0,1], [kf, r0, offz ])
            self.system.addForce(zrestrforce)

            
            #    atm_utils.addRestraintForce(lig_cm_particles = lig_atom_restr,
            #                            rcpt_cm_particles = rcpt_atom_restr,
            #                            kfcm = kf,
            #                            tolcm = r0,
            #                            offset = ligoffset)

        #indexes of the atoms whose position is restrained near the initial positions
        #by a flat-bottom harmonic potential. 
        posrestr_atoms_list = self.keywords.get('POS_RESTRAINED_ATOMS')
        if posrestr_atoms_list:
            posrestr_atoms = [int(i) for i in posrestr_atoms_list]
            fc = float(self.keywords.get('POSRE_FORCE_CONSTANT')) * kilocalories_per_mole
            tol = float(self.keywords.get('POSRE_TOLERANCE')) * angstrom
            atm_utils.addPosRestraints(posrestr_atoms, self.positions, fc, tol)
            
        #these define the state and will be overriden in set_state()
        temperature = 300 * kelvin
        lmbd = 0.0
        lambda1 = lmbd
        lambda2 = lmbd
        alpha = 0.0 / kilocalories_per_mole
        u0 = 0.0 * kilocalories_per_mole
        w0coeff = 0.0 * kilocalories_per_mole

        #soft-core parameters are fixed (the same in all states)
        umsc = float(self.keywords.get('UMAX')) * kilocalories_per_mole
        ubcore = self.keywords.get('UBCORE')
        if ubcore:
            ubcore = float(ubcore) * kilocalories_per_mole
        else:
            ubcore = 0.0 * kilocalories_per_mole
        acore = float(self.keywords.get('ACORE'))
        
        if not (self.keywords.get('DISPLACEMENT') is None):
            self.displ = [float(displ) for displ in self.keywords.get('DISPLACEMENT').split(',')]*angstrom
        else:
            msg = "Error: DISPLACEMENT is required"
            self._exit(msg)
        
        #create ATM Force
        self.atmforce = ATMMetaForce(lambda1, lambda2,  alpha * kilojoules_per_mole, u0/kilojoules_per_mole, w0coeff/kilojoules_per_mole, umsc/kilojoules_per_mole, ubcore/kilojoules_per_mole, acore)

        for i in range(self.topology.getNumAtoms()):
            self.atmforce.addParticle(i, 0., 0., 0.)
        for i in lig_atoms:
            self.atmforce.setParticleParameters(i, i, self.displ[0], self.displ[1], self.displ[2] )
        self.atmforce.setForceGroup(3)
        self.system.addForce(self.atmforce)
        
        #add barostat
        barostat = MonteCarloBarostat(1*bar, temperature)
        barostat.setForceGroup(1)                                                         
        barostat.setFrequency(0)#disabled
        self.system.addForce(barostat)

        #hack to store ASyncRE quantities in the openmm State
        sforce = mm.CustomBondForce("1")
        for name in self.parameter:
            sforce.addGlobalParameter(self.parameter[name], 0)
        sforce.setForceGroup(1)
        self.system.addForce(sforce)
        
        #temperature = int(self.keywords.get('TEMPERATURES')) * kelvin
        self.integrator = LangevinIntegrator(temperature/kelvin, self.frictionCoeff/(1/picosecond), self.MDstepsize/ picosecond )
        self.integrator.setIntegrationForceGroups({1,3})

        #these are the global parameters specified in the cntl files that need to be reset after reading the first
        #configuration
        self.cparams["ATMUmax"] = umsc/kilojoules_per_mole
        self.cparams["ATMUbcore"] = ubcore/kilojoules_per_mole
        self.cparams["ATMAcore"] = acore

class openmm_job_AmberABFE_zrestr(openmm_job_ATM):
    def __init__(self, command_file, options):
        super().__init__(command_file, options)

        prmtopfile = self.basename + ".prmtop"
        crdfile = self.basename + ".inpcrd"

        if not self.stateparams:
            self._buildStates()
        
        #builds service worker for replicas use
        service_ommsys = OMMSystemAmberABFE_zrestr(self.basename, self.keywords, prmtopfile, crdfile)
        self.service_worker = OMMWorkerATM(self.basename, service_ommsys, self.keywords, compute = False)
        #creates openmm replica objects
        self.openmm_replicas = []
        for i in range(self.nreplicas):
            replica = OMMReplicaATM(i, self.basename, self.service_worker)
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
            ommsys = OMMSystemAmberABFE_zrestr(self.basename, self.keywords, prmtopfile, crdfile)
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