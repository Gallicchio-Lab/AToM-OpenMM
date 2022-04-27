from __future__ import print_function
from __future__ import division
"""
Collects all of the ways that openmm systems are loaded
"""
import os, re, sys, time, shutil, copy, random, signal
import multiprocessing as mp
#from multiprocessing import Process, Queue, Event
import logging

from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from datetime import datetime
from configobj import ConfigObj

from atmmetaforce import *

class OMMSystem(object):
    def __init__(self, basename, keywords, logger):
        self.system = None
        self.topology = None
        self.positions = None
        self.boxvectors = None
        self.integrator = None
        self.keywords = keywords
        self.basename = basename
        self.logger = logger

        #parameters stored in the openmm state
        self.parameter = {}
        self.parameter['stateid'] = 'REStateId'
        self.parameter['cycle'] = 'RECycle'
        self.parameter['mdsteps'] = 'REMDSteps'
        #more ATM property names are in atmmetaforce

        #parameters from the cntl file
        self.cparams = {}

        self.frictionCoeff = float(self.keywords.get('FRICTION_COEFF')) / picosecond
        self.MDstepsize = float(self.keywords.get('TIME_STEP')) * picosecond

    def _exit(message):
        """Print and flush a message to stdout and then exit."""
        self.logger.error(message)
        sys.stdout.flush()
        sys.exit(1)

#Temperature RE
class OMMSystemAmberTRE(OMMSystem):
    def __init__(self, basename, keywords, prmtopfile, crdfile, logger):
        super().__init__(basename, keywords, logger)
        self.prmtopfile = prmtopfile
        self.crdfile = crdfile
        self.parameter['temperature'] = 'RETemperature'
        self.parameter['potential_energy'] = 'REPotEnergy'
        
    def create_system(self):

        self.prmtop = AmberPrmtopFile(self.prmtopfile)
        self.inpcrd = AmberInpcrdFile(self.crdfile)
        self.system = self.prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                                          constraints=HBonds)
        self.topology = self.prmtop.topology
        self.positions = self.inpcrd.positions

        #the temperature defines the state and will be overriden in set_state()
        temperature = 300 * kelvin

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
        
        self.integrator = LangevinIntegrator(temperature/kelvin, self.frictionCoeff/(1/picosecond), self.MDstepsize/ picosecond )

class OMMSystemAmberABFE(OMMSystem):
    def __init__(self, basename, keywords, prmtopfile, crdfile, logger):
        super().__init__(basename, keywords, logger)
        self.prmtopfile = prmtopfile
        self.crdfile = crdfile
        self.parameter['temperature'] = 'RETemperature'
        self.parameter['potential_energy'] = 'REPotEnergy'
        self.parameter['perturbation_energy'] = 'REPertEnergy'
        self.parameter['atmintermediate'] = 'REAlchemicalIntermediate'
        self.atmforce = None
        self.lig_atoms = None
        self.displ = None
        
    def create_system(self):

        self.prmtop = AmberPrmtopFile(self.prmtopfile)
        self.inpcrd = AmberInpcrdFile(self.crdfile)
        self.system = self.prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                                          constraints=HBonds)
        self.topology = self.prmtop.topology
        self.positions = self.inpcrd.positions
        self.boxvectors = self.inpcrd.boxVectors
        
        atm_utils = ATMMetaForceUtils(self.system)
        
        lig_atoms_in = self.keywords.get('LIGAND_ATOMS')   #indexes of ligand atoms
        if lig_atoms_in is not None:
            self.lig_atoms = [int(i) for i in lig_atoms_in]
        else:
            msg = "Error: LIGAND_ATOMS is required"
            self._exit(msg)


        #CM-CM Vsite restraints
        cm_lig_atoms = self.keywords.get('LIGAND_CM_ATOMS')   #indexes of ligand atoms for CM-CM Vsite restraint
        if cm_lig_atoms is not None:
            lig_atom_restr = [int(i) for i in cm_lig_atoms]
        else:
            lig_atom_restr = None
        cm_rcpt_atoms = self.keywords.get('RCPT_CM_ATOMS')   #indexes of rcpt atoms for CM-CM Vsite restraint
        if cm_rcpt_atoms is not None:
            rcpt_atom_restr = [int(i) for i in cm_rcpt_atoms]
        else:
            rcpt_atom_restr = None
        cmrestraints_present = (cm_rcpt_atoms is not None) and (cm_lig_atoms is not None)
        self.vsiterestraintForce = None
        if cmrestraints_present:
            cmkf = float(self.keywords.get('CM_KF'))
            kf = cmkf * kilocalorie_per_mole/angstrom**2 #force constant for Vsite CM-CM restraint
            cmtol = float(self.keywords.get('CM_TOL'))
            r0 = cmtol * angstrom #radius of Vsite sphere
            ligoffset = self.keywords.get('LIGOFFSET')
            if ligoffset is not None:
                ligoffset = [float(offset) for offset in ligoffset.split(',')]*angstrom
            self.vsiterestraintForce = atm_utils.addRestraintForce(lig_cm_particles = lig_atom_restr,
                                                                   rcpt_cm_particles = rcpt_atom_restr,
                                                                   kfcm = kf,
                                                                   tolcm = r0,
                                                                   offset = ligoffset)

        #orientational VSite restraints
        #the indexes of the groups of atoms that define the internal reference frame of the ligand
        lig_frame_groups = None
        lig_frame_groups_inp = self.keywords.get('LIGAND_VSITE_FRAMEGROUPS')
        if lig_frame_groups_inp is not None:
            lig_frame_groups = []
            for i in range(3):
                lig_frame_groups.append([int(j) for j in lig_frame_groups_inp[str(i)]])
        #the indexes of the groups of atoms that define the internal reference frame of the receptor
        rcpt_frame_groups = None
        rcpt_frame_groups_inp = self.keywords.get('RCPT_VSITE_FRAMEGROUPS')
        if rcpt_frame_groups_inp is not None:
            rcpt_frame_groups = []
            for i in range(3):
                rcpt_frame_groups.append([int(j) for j in rcpt_frame_groups_inp[str(i)]])
        if (lig_frame_groups is not None) and (rcpt_frame_groups is not None):
             kftheta = self.keywords.get('VSITE_KFTHETA')
             theta0 = self.keywords.get('VSITE_THETA0')
             thetatol = self.keywords.get('VSITE_THETATOL')
             kfphi = self.keywords.get('VSITE_KFPHI')
             phi0 = self.keywords.get('VSITE_PHI0')
             phitol = self.keywords.get('VSITE_PHITOL')
             kfpsi = self.keywords.get('VSITE_KFPSI')
             psi0 = self.keywords.get('VSITE_PSI0')
             psitol = self.keywords.get('VSITE_PSITOL')
             kftheta = kftheta if kftheta is None else float(kftheta)*kilocalories_per_mole
             theta0 = theta0 if theta0 is None else float(theta0)*degrees
             thetatol = thetatol if thetatol is None else float(thetatol)*degrees
             kfphi = kfphi if kfphi is None else float(kfphi)*(kilocalories_per_mole/degrees**2)
             phi0 = phi0 if phi0 is None else float(phi0)*degrees
             phitol = phitol if phitol is None else float(phitol)*degrees
             kfpsi = kfpsi if kfpsi is None else float(kfpsi)*(kilocalories_per_mole/degrees**2)
             psi0 = psi0 if psi0 is None else float(psi0)*degrees
             psitol = psitol if psitol is None else float(psitol)*degrees
             atm_utils.addVsiteRestraintForceCMAngles(lig_frame_groups, rcpt_frame_groups, 
                                                      kftheta, theta0, thetatol,
                                                      kfphi, phi0, phitol,
                                                      kfpsi, psi0, psitol)

        #indexes of the atoms whose position is restrained near the initial positions
        #by a flat-bottom harmonic potential. 
        posrestr_atoms_list = self.keywords.get('POS_RESTRAINED_ATOMS')
        self.posrestrForce = None
        if posrestr_atoms_list is not None:
            posrestr_atoms = [int(i) for i in posrestr_atoms_list]
            fc = float(self.keywords.get('POSRE_FORCE_CONSTANT')) * kilocalorie_per_mole
            tol = float(self.keywords.get('POSRE_TOLERANCE')) * angstrom
            self.posrestrForce = atm_utils.addPosRestraints(posrestr_atoms, self.positions, fc, tol)
            
        #these define the state and will be overriden in set_state()
        temperature = 300 * kelvin
        lmbd = 0.0
        lambda1 = lmbd
        lambda2 = lmbd
        alpha = 0.0 / kilocalorie_per_mole
        u0 = 0.0 * kilocalorie_per_mole
        w0coeff = 0.0 * kilocalorie_per_mole
        alchemical_direction = 1.0

        #soft-core parameters are fixed (the same in all states)
        umsc = float(self.keywords.get('UMAX')) * kilocalorie_per_mole
        ubcore = self.keywords.get('UBCORE')
        if ubcore is not None:
            ubcore = float(ubcore) * kilocalorie_per_mole
        else:
            ubcore = 0.0 * kilocalorie_per_mole
        acore = float(self.keywords.get('ACORE'))

        if not (self.keywords.get('DISPLACEMENT') is None):
            self.displ = [float(displ) for displ in self.keywords.get('DISPLACEMENT').split(',')]*angstrom
        else:
            msg = "Error: DISPLACEMENT is required"
            self._exit(msg)
        
        #create ATM Force
        self.atmforce = ATMMetaForce(lambda1, lambda2,  alpha * kilojoules_per_mole, u0/kilojoules_per_mole, w0coeff/kilojoules_per_mole, umsc/kilojoules_per_mole, ubcore/kilojoules_per_mole, acore, alchemical_direction)

        for i in range(self.topology.getNumAtoms()):
            self.atmforce.addParticle(i, 0., 0., 0.)
        for i in self.lig_atoms:
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

        #these are the global parameters specified in the cntl files that need to be reset
        #by the worker after reading the first configuration
        self.cparams["ATMUmax"] = umsc/kilojoules_per_mole
        self.cparams["ATMUbcore"] = ubcore/kilojoules_per_mole
        self.cparams["ATMAcore"] = acore

        
class OMMSystemAmberRBFE(OMMSystem):
    def __init__(self, basename, keywords, prmtopfile, crdfile, logger):
        super().__init__(basename, keywords, logger)
        self.prmtopfile = prmtopfile
        self.crdfile = crdfile
        self.parameter['temperature'] = 'RETemperature'
        self.parameter['potential_energy'] = 'REPotEnergy'
        self.parameter['perturbation_energy'] = 'REPertEnergy'
        self.parameter['atmintermediate'] = 'REAlchemicalIntermediate'
        self.atmforce = None
        self.lig1_atoms = None
        self.lig2_atoms = None
        self.displ = None
        
    def create_system(self):

        self.prmtop = AmberPrmtopFile(self.prmtopfile)
        self.inpcrd = AmberInpcrdFile(self.crdfile)
        self.system = self.prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                                          constraints=HBonds)
        self.topology = self.prmtop.topology
        self.positions = self.inpcrd.positions

        atm_utils = ATMMetaForceUtils(self.system)

        lig1_atoms_in = self.keywords.get('LIGAND1_ATOMS')   #indexes of ligand1 atoms
        lig2_atoms_in = self.keywords.get('LIGAND2_ATOMS')   #indexes of ligand2 atoms
        if lig1_atoms_in is not None:
            self.lig1_atoms = [int(i) for i in lig1_atoms_in]
        else:
            msg = "Error: LIGAND1_ATOMS is required"
            self._exit(msg)        
        if lig2_atoms_in is not None:
            self.lig2_atoms = [int(i) for i in lig2_atoms_in ]
        else:
            msg = "Error: LIGAND2_ATOMS is required"
            self._exit(msg)
        
        #ligand 1 Vsite restraint
        cm_lig1_atoms = self.keywords.get('LIGAND1_CM_ATOMS')   #indexes of ligand atoms for CM-CM Vsite restraint
        if cm_lig1_atoms is not None:
            lig1_atom_restr = [int(i) for i in cm_lig1_atoms]
        else:
            lig1_atom_restr = None
        
        #ligand 2 Vsite restraint
        cm_lig2_atoms = self.keywords.get('LIGAND2_CM_ATOMS')   #indexes of ligand atoms for CM-CM Vsite restraint
        if cm_lig2_atoms is not None:
            lig2_atom_restr = [int(i) for i in cm_lig2_atoms]
        else:
            lig2_atom_restr = None
        
        #Vsite restraint receptor atoms
        cm_rcpt_atoms = self.keywords.get('RCPT_CM_ATOMS')   #indexes of rcpt atoms for CM-CM Vsite restraint
        if cm_rcpt_atoms is not None:
            rcpt_atom_restr = [int(i) for i in cm_rcpt_atoms]
        else:
            rcpt_atom_restr = None

        #set displacements and offsets for ligand 1 and ligand 2
        if self.keywords.get('DISPLACEMENT') is not None:
            self.displ = [float(displ) for displ in self.keywords.get('DISPLACEMENT').split(',')]*angstrom
            self.lig1offset = [float(0.0*offset) for offset in self.displ/angstrom]*angstrom
            self.lig2offset = [float(offset) for offset in self.displ/angstrom]*angstrom
        else:
            msg = "DISPLACEMENT is required"
            self._exit(msg)
            
        cmrestraints_present = (rcpt_atom_restr is not None) and (lig1_atom_restr is not None) and (lig2_atom_restr is not None)

        self.vsiterestraintForce1 = None
        self.vsiterestraintForce2 = None
        if cmrestraints_present:
            cmkf = float(self.keywords.get('CM_KF'))
            kf = cmkf * kilocalorie_per_mole/angstrom**2 #force constant for Vsite CM-CM restraint
            cmtol = float(self.keywords.get('CM_TOL'))
            r0 = cmtol * angstrom #radius of Vsite sphere            

            #Vsite restraints for ligands 1 and 2
            self.vsiterestraintForce1 = atm_utils.addRestraintForce(lig_cm_particles = lig1_atom_restr,
                                        rcpt_cm_particles = rcpt_atom_restr,
                                        kfcm = kf,
                                        tolcm = r0,
                                        offset = self.lig1offset)
            self.vsiterestraintForce2 = atm_utils.addRestraintForce(lig_cm_particles = lig2_atom_restr,
                                        rcpt_cm_particles = rcpt_atom_restr,
                                        kfcm = kf,
                                        tolcm = r0,
                                        offset = self.lig2offset)

        #orientational VSite restraints
        #the indexes of the groups of atoms that define the internal reference frame of the ligand
        lig1_frame_groups = None
        lig1_frame_groups_inp = self.keywords.get('LIGAND1_VSITE_FRAMEGROUPS')
        if lig1_frame_groups_inp is not None:
            lig1_frame_groups = []
            for i in range(3):
                lig1_frame_groups.append([int(j) for j in lig1_frame_groups_inp[str(i)]])
        lig2_frame_groups = None
        lig2_frame_groups_inp = self.keywords.get('LIGAND2_VSITE_FRAMEGROUPS')
        if lig2_frame_groups_inp is not None:
            lig2_frame_groups = []
            for i in range(3):
                lig2_frame_groups.append([int(j) for j in lig2_frame_groups_inp[str(i)]])
        #the indexes of the groups of atoms that define the internal reference frame of the receptor
        rcpt_frame_groups = None
        rcpt_frame_groups_inp = self.keywords.get('RCPT_VSITE_FRAMEGROUPS')
        if rcpt_frame_groups_inp is not None:
            rcpt_frame_groups = []
            for i in range(3):
                rcpt_frame_groups.append([int(j) for j in rcpt_frame_groups_inp[str(i)]])
        if rcpt_frame_groups is not None:
            kftheta = self.keywords.get('VSITE_KFTHETA_LIG1')
            theta0 = self.keywords.get('VSITE_THETA0_LIG1')
            thetatol = self.keywords.get('VSITE_THETATOL_LIG1')
            kfphi = self.keywords.get('VSITE_KFPHI_LIG1')
            phi0 = self.keywords.get('VSITE_PHI0_LIG1')
            phitol = self.keywords.get('VSITE_PHITOL_LIG1')
            kfpsi = self.keywords.get('VSITE_KFPSI_LIG1')
            psi0 = self.keywords.get('VSITE_PSI0_LIG1')
            psitol = self.keywords.get('VSITE_PSITOL_LIG1')
            kftheta = kftheta if kftheta is None else float(kftheta)*kilocalories_per_mole
            theta0 = theta0 if theta0 is None else float(theta0)*degrees
            thetatol = thetatol if thetatol is None else float(thetatol)*degrees
            kfphi = kfphi if kfphi is None else float(kfphi)*(kilocalories_per_mole/degrees**2)
            phi0 = phi0 if phi0 is None else float(phi0)*degrees
            phitol = phitol if phitol is None else float(phitol)*degrees
            kfpsi = kfpsi if kfpsi is None else float(kfpsi)*(kilocalories_per_mole/degrees**2)
            psi0 = psi0 if psi0 is None else float(psi0)*degrees
            psitol = psitol if psitol is None else float(psitol)*degrees
            if lig1_frame_groups is not None:
                atm_utils.addVsiteRestraintForceCMAngles(lig1_frame_groups, rcpt_frame_groups,
                                                         kftheta, theta0, thetatol,
                                                         kfphi, phi0, phitol,
                                                         kfpsi, psi0, psitol)
            kftheta = self.keywords.get('VSITE_KFTHETA_LIG2')
            theta0 = self.keywords.get('VSITE_THETA0_LIG2')
            thetatol = self.keywords.get('VSITE_THETATOL_LIG2')
            kfphi = self.keywords.get('VSITE_KFPHI_LIG2')
            phi0 = self.keywords.get('VSITE_PHI0_LIG2')
            phitol = self.keywords.get('VSITE_PHITOL_LIG2')
            kfpsi = self.keywords.get('VSITE_KFPSI_LIG2')
            psi0 = self.keywords.get('VSITE_PSI0_LIG2')
            psitol = self.keywords.get('VSITE_PSITOL_LIG2')
            kftheta = kftheta if kftheta is None else float(kftheta)*kilocalories_per_mole
            theta0 = theta0 if theta0 is None else float(theta0)*degrees
            thetatol = thetatol if thetatol is None else float(thetatol)*degrees
            kfphi = kfphi if kfphi is None else float(kfphi)*(kilocalories_per_mole/degrees**2)
            phi0 = phi0 if phi0 is None else float(phi0)*degrees
            phitol = phitol if phitol is None else float(phitol)*degrees
            kfpsi = kfpsi if kfpsi is None else float(kfpsi)*(kilocalories_per_mole/degrees**2)
            psi0 = psi0 if psi0 is None else float(psi0)*degrees
            psitol = psitol if psitol is None else float(psitol)*degrees
            if lig2_frame_groups is not None:
                atm_utils.addVsiteRestraintForceCMAngles(lig2_frame_groups, rcpt_frame_groups,
                                                         kftheta, theta0, thetatol,
                                                         kfphi, phi0, phitol,
                                                         kfpsi, psi0, psitol)

        #reference atoms for alignment force
        refatoms1_cntl = self.keywords.get('ALIGN_LIGAND1_REF_ATOMS')
        self.refatoms1 = [int(refatoms1) for refatoms1 in refatoms1_cntl]
        lig1_ref_atoms  = [ self.refatoms1[i]+self.lig1_atoms[0] for i in range(3)]
        
        refatoms2_cntl = self.keywords.get('ALIGN_LIGAND2_REF_ATOMS')
        self.refatoms2 = [int(refatoms2) for refatoms2 in refatoms2_cntl]
        lig2_ref_atoms  = [ self.refatoms2[i]+self.lig2_atoms[0] for i in range(3)]
        
        #add alignment force
        atm_utils.addAlignmentForce(liga_ref_particles = lig1_ref_atoms,
                                    ligb_ref_particles = lig2_ref_atoms,
                                    kfdispl = float(self.keywords.get('ALIGN_KF_SEP'))*kilocalorie_per_mole/angstrom**2,
                                    ktheta = float(self.keywords.get('ALIGN_K_THETA'))*kilocalorie_per_mole,
                                    kpsi = float(self.keywords.get('ALIGN_K_PSI'))*kilocalorie_per_mole,
                                    offset = self.lig2offset)

        #indexes of the atoms whose position is restrained near the initial positions
        #by a flat-bottom harmonic potential. 
        posrestr_atoms_list = self.keywords.get('POS_RESTRAINED_ATOMS')
        self.posrestrForce = None
        if posrestr_atoms_list is not None:
            posrestr_atoms = [int(i) for i in posrestr_atoms_list]
            fc = float(self.keywords.get('POSRE_FORCE_CONSTANT')) * kilocalorie_per_mole
            tol = float(self.keywords.get('POSRE_TOLERANCE')) * angstrom
            self.posrestrForce = atm_utils.addPosRestraints(posrestr_atoms, self.positions, fc, tol)
            
        #these define the state and will be overriden in set_state()
        temperature = 300 * kelvin
        lmbd = 0.0
        lambda1 = lmbd
        lambda2 = lmbd
        alpha = 0.0 / kilocalorie_per_mole
        u0 = 0.0 * kilocalorie_per_mole
        w0coeff = 0.0 * kilocalorie_per_mole

        #soft-core parameters are fixed (the same in all states)
        umsc = float(self.keywords.get('UMAX')) * kilocalorie_per_mole
        ubcore = self.keywords.get('UBCORE')
        if ubcore:
            ubcore = float(ubcore) * kilocalorie_per_mole
        else:
            ubcore = 0.0 * kilocalorie_per_mole
        acore = float(self.keywords.get('ACORE'))

        #create ATM Force
        self.atmforce = ATMMetaForce(lambda1, lambda2,  alpha * kilojoules_per_mole, u0/kilojoules_per_mole, w0coeff/kilojoules_per_mole, umsc/kilojoules_per_mole, ubcore/kilojoules_per_mole, acore)

        for i in range(self.topology.getNumAtoms()):
            self.atmforce.addParticle(i, 0., 0., 0.)
        for i in self.lig1_atoms:
            self.atmforce.setParticleParameters(i, i, self.displ[0], self.displ[1], self.displ[2] )
        for i in self.lig2_atoms:
            self.atmforce.setParticleParameters(i, i, -self.displ[0], -self.displ[1], -self.displ[2] )

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
