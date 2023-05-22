from __future__ import print_function
from __future__ import division
"""
Collects all of the ways that openmm systems are loaded
"""
import os, re, sys, time, shutil, copy, random, signal
import numpy as np
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

# OpenMM's MTSLangevinIntegrator does not have a setTemperature method
class ATMMTSLangevinIntegrator(MTSLangevinIntegrator):
    def setTemperature(self, temperature):
        self.setGlobalVariableByName('kT', MOLAR_GAS_CONSTANT_R*temperature)

class OMMSystem(object):
    def __init__(self, basename, keywords, logger):
        self.system = None
        self.topology = None
        self.positions = None
        self.boxvectors = None
        self.integrator = None
        self.barostat = None
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

        self.atmforcegroup = 2
        self.nonbondedforcegroup = 1
        self.metaDforcegroup = 3

        self.doMetaD = False

    def _exit(message):
        """Print and flush a message to stdout and then exit."""
        self.logger.error(message)
        sys.stdout.flush()
        sys.exit(1)

class OMMSystemAmber(OMMSystem):
    def __init__(self, basename, keywords, prmtopfile, crdfile, logger):
        super().__init__(basename, keywords, logger)
        self.prmtopfile = prmtopfile
        self.crdfile = crdfile
        self.parameter['temperature'] = 'RETemperature'
        self.parameter['potential_energy'] = 'REPotEnergy'

    def load_amber_system(self):
        """
        sets the value of
        prmtop : Amber topology file
        inpcrd : Amber coordinates
        system : creates a OpenMM system with topology, coordinates
        topology : defines the OpenMM topology of the system
        positions : defines positions of all atoms of the system in OpenMM
        boxvectors : stores the dimension of the simulation box

        """
        self.prmtop = AmberPrmtopFile(self.prmtopfile)
        self.inpcrd = AmberInpcrdFile(self.crdfile)
        if self.keywords.get('HMASS') is not None:
            hmass = float(self.keywords.get('HMASS'))*amu
        else:
            hmass = 1.0*amu
        self.system = self.prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=0.9*nanometer,
                                               constraints=HBonds, hydrogenMass = hmass)
        self.topology = self.prmtop.topology
        self.positions = self.inpcrd.positions
        self.boxvectors = self.inpcrd.boxVectors

        if nnp_model := self.keywords.get('NNP_MODEL'):
            from openmmml import MLPotential
            self.logger.info(f'NNP model: {nnp_model}')
            nnp = MLPotential(nnp_model)
            nnp_atoms = list(map(int, self.keywords['LIGAND1_ATOMS'] + self.keywords['LIGAND2_ATOMS']))
            self.logger.info(f'NNP atoms: {nnp_atoms}')
            self.logger.info(f'Number of NNP atoms: {len(nnp_atoms)}')
            assert len(nnp_atoms) > 0
            self.system = nnp.createMixedSystem(self.topology, self.system, nnp_atoms, removeConstraints=False)

    def set_barostat(self,temperature,pressure,frequency):
        """
        sets the system Barostat; Currently applies the MonteCarlo Barostat

        Requires: self,
        temperature
        pressure : eg. 1*bar
        frequency : 0 - disable the barostat

        """
        self.barostat = MonteCarloBarostat(pressure, temperature)
        self.barostat.setFrequency(frequency)
        self.system.addForce(self.barostat)

    def set_integrator(self, temperature, frictionCoeff, MDstepsize, defaultMDstepsize = 0.001*picosecond):
        #place non-bonded force in group 1, assume all other bonded forces are in group 0
        nonbonded = [f for f in self.system.getForces() if isinstance(f, NonbondedForce)][0]
        nonbonded.setForceGroup(self.nonbondedforcegroup)
        #set the multiplicity of the calculation of bonded forces so that they are evaluated at least once every 1 fs (default time-step)
        bonded_frequency = max(1, int(round(MDstepsize/defaultMDstepsize)))
        self.logger.info("Running with a %f fs time-step with bonded forces integrated %d times per time-step" % (MDstepsize/femtosecond, bonded_frequency))
        if doMetaD:
            fgroups = [(0,bonded_frequency), (self.metaDforcegroup, bonded_frequency), (self.nonbondedforcegroup,1)]
        else:
            fgroups = [(0,bonded_frequency), (self.nonbondedforcegroup,1)]
        self.integrator = ATMMTSLangevinIntegrator(temperature, frictionCoeff, MDstepsize, fgroups)
        self.integrator.setConstraintTolerance(0.00001)

    def set_positional_restraints(self):
        #indexes of the atoms whose position is restrained near the initial positions
        #by a flat-bottom harmonic potential. 
        posrestr_atoms_list = self.keywords.get('POS_RESTRAINED_ATOMS')
        self.posrestrForce = None
        if posrestr_atoms_list is not None:
            posrestr_atoms = [int(i) for i in posrestr_atoms_list]
            fc = float(self.keywords.get('POSRE_FORCE_CONSTANT')) * (kilocalorie_per_mole/angstrom**2)
            tol = float(self.keywords.get('POSRE_TOLERANCE')) * angstrom
            self.posrestrForce = self.atm_utils.addPosRestraints(posrestr_atoms, self.positions, fc, tol)

    def set_torsion_metaDbias(self,temperature):
        if self.keywords.get('METADBIAS_DIR') == None:
            return
        bias_dirs = self.keywords.get('METADBIAS_DIR')
        bias_offsets = self.keywords.get('METADBIAS_IDXOFFSET')

        for mdir,offset in zip(bias_dirs,bias_offsets) :
            cntlfile = "%s/%s.cntl" % (mdir, mdir)
            keywords  = ConfigObj(cntlfile, file_error = True)

            #metadynamics settings
            bias_factor = float(keywords.get('METADBIAS_FACTOR')) # this is (T+DeltaT)/T
            bias_height = float(keywords.get('METADBIAS_GHEIGHT')) * kilocalorie_per_mole #height of each gaussian
            bias_frequency = int(keywords.get('METADBIAS_FREQUENCY')) #steps in between gaussian depositions
            bias_savefrequency = int(keywords.get('METADBIAS_SAVEFREQUENCY')) #steps in between checkpointing of bias potential

            #bias force settings
            torsions = keywords.get('METADBIAS_TORSIONS')
            ndim = len(torsions.keys())

            gaussian_width = keywords.get('METADBIAS_GWIDTH')
            angle_min = keywords.get('METADBIAS_MINANGLE')
            angle_max = keywords.get('METADBIAS_MAXANGLE')
            ngrid = keywords.get('METADBIAS_NGRID')
            periodic = keywords.get('METADBIAS_PERIODIC')

            torForce = []
            biasvar = []

            for t in range(ndim):
                torForce.append(mm.CustomTorsionForce("theta"))
                p = torsions[str(t)]
                gw = float(gaussian_width[t])*kilocalorie_per_mole
                amin = float(angle_min[t])*degrees
                amax = float(angle_max[t])*degrees
                per = int(periodic[t]) > 0
                ng = int(ngrid[t])
                dp = int(offset)
                torForce[t].addTorsion(int(p[0])+dp, int(p[1])+dp, int(p[2])+dp, int(p[3])+dp)
                biasvar.append(BiasVariable(torForce[t], amin, amax, gw, per, ng))

            metaD = Metadynamics(self.system, biasvar, temperature, bias_factor, bias_height, bias_frequency, bias_savefrequency, mdir)
            metaD._force.setForceGroup(self.metaDforcegroup)
        self.doMetaD = True

#Temperature RE
class OMMSystemAmberTRE(OMMSystemAmber):
    def create_system(self):
        self.load_amber_system()
        #the temperature defines the state and will be overriden in set_state()
        temperature = 300 * kelvin
        #set barostat
        self.set_barostat(temperature,1*bar,900000000)

        #hack to store ASyncRE quantities in the openmm State
        sforce = mm.CustomBondForce("1")
        for name in self.parameter:
            sforce.addGlobalParameter(self.parameter[name], 0)
        self.system.addForce(sforce)

        self.set_integrator(temperature, self.frictionCoeff, self.MDstepsize)

class OMMSystemAmberABFE(OMMSystemAmber):
    def __init__(self, basename, keywords, prmtopfile, crdfile, logger):
        super().__init__(basename, keywords, prmtopfile, crdfile, logger)

        self.parameter['perturbation_energy'] = 'REPertEnergy'
        self.parameter['atmintermediate'] = 'REAlchemicalIntermediate'
        self.parameter['bias_energy'] = 'BiasEnergy'
        self.atmforce = None
        self.lig_atoms = None
        self.displ = None

    def set_ligand_atoms(self):
        lig_atoms_in = self.keywords.get('LIGAND_ATOMS')   #indexes of ligand atoms
        if lig_atoms_in is not None:
            self.lig_atoms = [int(i) for i in lig_atoms_in]
        else:
            msg = "Error: LIGAND_ATOMS is required"
            self._exit(msg)

    def set_vsite_restraints(self):
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
            self.vsiterestraintForce = self.atm_utils.addRestraintForce(lig_cm_particles = lig_atom_restr,
                                                                   rcpt_cm_particles = rcpt_atom_restr,
                                                                   kfcm = kf,
                                                                   tolcm = r0,
                                                                   offset = ligoffset)

    def set_orientation_restraints(self):
        #orientation VSite restraints
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
             self.atm_utils.addVsiteRestraintForceCMAngles(lig_frame_groups, rcpt_frame_groups, 
                                                      kftheta, theta0, thetatol,
                                                      kfphi, phi0, phitol,
                                                      kfpsi, psi0, psitol)

    def set_integrator(self, temperature, frictionCoeff, MDstepsize, defaultMDstepsize = 0.001*picosecond):
        #set the multiplicity of the calculation of bonded forces so that they are evaluated at least once every 1 fs (default time-step)
        bonded_frequency = max(1, int(round(MDstepsize/defaultMDstepsize)))
        self.logger.info("Running with a %f fs time-step with bonded forces integrated %d times per time-step" % (MDstepsize/femtosecond, bonded_frequency))
        if self.doMetaD:
            fgroups = [(0,bonded_frequency), (self.metaDforcegroup, bonded_frequency), (self.atmforcegroup,1)]
        else:
            fgroups = [(0,bonded_frequency), (self.atmforcegroup,1)]
        self.integrator = ATMMTSLangevinIntegrator(temperature, frictionCoeff, MDstepsize, fgroups )
        self.integrator.setConstraintTolerance(0.00001)

    def set_atmforce(self):
        #these define the state and will be overriden in set_state()
        lmbd = 0.0
        lambda1 = lmbd
        lambda2 = lmbd
        alpha = 0.0 / kilocalorie_per_mole
        u0 = 0.0 * kilocalorie_per_mole
        w0coeff = 0.0 * kilocalorie_per_mole
        direction = 1.0

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
        self.atm_utils.setNonbondedForceGroup(self.nonbondedforcegroup)
        atmvariableforcegroups = [self.nonbondedforcegroup]
        self.atmforce = ATMMetaForce(lambda1, lambda2,  alpha * kilojoules_per_mole, u0/kilojoules_per_mole, w0coeff/kilojoules_per_mole, umsc/kilojoules_per_mole, ubcore/kilojoules_per_mole, acore, direction, atmvariableforcegroups )

        for i in range(self.topology.getNumAtoms()):
            self.atmforce.addParticle(i, 0., 0., 0.)
        for i in self.lig_atoms:
            self.atmforce.setParticleParameters(i, i, self.displ[0], self.displ[1], self.displ[2] )
        self.atmforce.setForceGroup(self.atmforcegroup)
        self.system.addForce(self.atmforce)
        #these are the global parameters specified in the cntl files that need to be reset
        #by the worker after reading the first configuration
        self.cparams["ATMUmax"] = umsc/kilojoules_per_mole
        self.cparams["ATMUbcore"] = ubcore/kilojoules_per_mole
        self.cparams["ATMAcore"] = acore

    def create_system(self):
        self.load_amber_system()

        self.atm_utils = ATMMetaForceUtils(self.system)

        self.set_ligand_atoms()

        self.set_vsite_restraints()

        self.set_orientation_restraints()

        self.set_positional_restraints()

        self.set_atmforce()

        #temperature is part of the state and is maybe overriden in set_state()
        temperature = 300 * kelvin

        #add barostat
        pressure=1*bar
        self.set_barostat(temperature,pressure,900000000)

        #hack to store ASyncRE quantities in the openmm State
        sforce = mm.CustomBondForce("1")
        for name in self.parameter:
            sforce.addGlobalParameter(self.parameter[name], 0)
        self.system.addForce(sforce)

        self.set_integrator(temperature, self.frictionCoeff, self.MDstepsize)


class OMMSystemAmberRBFE(OMMSystemAmber):
    def __init__(self, basename, keywords, prmtopfile, crdfile, logger):
        super().__init__(basename, keywords, prmtopfile, crdfile, logger)

        self.parameter['perturbation_energy'] = 'REPertEnergy'
        self.parameter['atmintermediate'] = 'REAlchemicalIntermediate'
        self.parameter['bias_energy'] = 'BiasEnergy'
        self.atmforce = None
        self.lig1_atoms = None
        self.lig2_atoms = None
        self.displ = None

    def set_ligand_atoms(self):
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

    def set_displacement(self):
        #set displacements and offsets for ligand 1 and ligand 2
        if self.keywords.get('DISPLACEMENT') is not None:
            self.displ = [float(displ) for displ in self.keywords.get('DISPLACEMENT').split(',')]*angstrom
            self.lig1offset = [float(0.0*offset) for offset in self.displ/angstrom]*angstrom
            self.lig2offset = [float(offset) for offset in self.displ/angstrom]*angstrom
        else:
            msg = "DISPLACEMENT is required"
            self._exit(msg)

    def set_vsite_restraints(self):
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
        if cm_rcpt_atoms is None:
            cm_rcpt_atoms = self.keywords.get('REST_LIGAND_CMREC_ATOMS')
        if cm_rcpt_atoms is not None:
            rcpt_atom_restr = [int(i) for i in cm_rcpt_atoms]
        else:
            rcpt_atom_restr = None

        cmrestraints_present = (rcpt_atom_restr is not None) and (lig1_atom_restr is not None) and (lig2_atom_restr is not None)

        self.vsiterestraintForce1 = None
        self.vsiterestraintForce2 = None
        if cmrestraints_present:
            cmkf = float(self.keywords.get('CM_KF'))
            kf = cmkf * kilocalorie_per_mole/angstrom**2 #force constant for Vsite CM-CM restraint
            cmtol = float(self.keywords.get('CM_TOL'))
            r0 = cmtol * angstrom #radius of Vsite sphere

            #Vsite restraints for ligands 1 and 2
            self.vsiterestraintForce1 = self.atm_utils.addRestraintForce(lig_cm_particles = lig1_atom_restr,
                                        rcpt_cm_particles = rcpt_atom_restr,
                                        kfcm = kf,
                                        tolcm = r0,
                                        offset = self.lig1offset)
            self.vsiterestraintForce2 = self.atm_utils.addRestraintForce(lig_cm_particles = lig2_atom_restr,
                                        rcpt_cm_particles = rcpt_atom_restr,
                                        kfcm = kf,
                                        tolcm = r0,
                                        offset = self.lig2offset)

    def set_orientation_restraints(self):
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
                self.atm_utils.addVsiteRestraintForceCMAngles(lig1_frame_groups, rcpt_frame_groups,
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
                self.atm_utils.addVsiteRestraintForceCMAngles(lig2_frame_groups, rcpt_frame_groups,
                                                         kftheta, theta0, thetatol,
                                                         kfphi, phi0, phitol,
                                                         kfpsi, psi0, psitol)

    def set_alignmentForce(self):
        """
        set reference atoms for adding the alignment force

        """
        refatoms1_cntl = self.keywords.get('ALIGN_LIGAND1_REF_ATOMS')
        refatoms2_cntl = self.keywords.get('ALIGN_LIGAND2_REF_ATOMS')

        if refatoms1_cntl == None or refatoms2_cntl == None:
            return

        self.refatoms1 = [int(refatoms1) for refatoms1 in refatoms1_cntl]
        lig1_ref_atoms  = [ self.refatoms1[i]+self.lig1_atoms[0] for i in range(3)]
        self.refatoms2 = [int(refatoms2) for refatoms2 in refatoms2_cntl]
        lig2_ref_atoms  = [ self.refatoms2[i]+self.lig2_atoms[0] for i in range(3)]

        #add alignment force
        self.atm_utils.addAlignmentForce(liga_ref_particles = lig1_ref_atoms,
                                    ligb_ref_particles = lig2_ref_atoms,
                                    kfdispl = float(self.keywords.get('ALIGN_KF_SEP'))*kilocalorie_per_mole/angstrom**2,
                                    ktheta = float(self.keywords.get('ALIGN_K_THETA'))*kilocalorie_per_mole,
                                    kpsi = float(self.keywords.get('ALIGN_K_PSI'))*kilocalorie_per_mole,
                                    offset = self.lig2offset)

    def set_integrator(self, temperature, frictionCoeff, MDstepsize, defaultMDstepsize = 0.001*picosecond):
        #set the multiplicity of the calculation of bonded forces so that they are evaluated at least once every 1 fs (default time-step)
        bonded_frequency = max(1, int(round(MDstepsize/defaultMDstepsize)))
        self.logger.info("Running with a %f fs time-step with bonded forces integrated %d times per time-step" % (MDstepsize/femtosecond, bonded_frequency))
        if self.doMetaD:
            fgroups = [(0,bonded_frequency), (self.metaDforcegroup, bonded_frequency), (self.atmforcegroup,1)]
        else:
            fgroups = [(0,bonded_frequency), (self.atmforcegroup,1)]
        self.integrator = ATMMTSLangevinIntegrator(temperature, frictionCoeff, MDstepsize, fgroups )
        self.integrator.setConstraintTolerance(0.00001)

    def set_atmforce(self):
        #these define the state and will be overriden in set_state()
        temperature = 300 * kelvin
        lmbd = 0.0
        lambda1 = lmbd
        lambda2 = lmbd
        alpha = 0.0 / kilocalorie_per_mole
        u0 = 0.0 * kilocalorie_per_mole
        w0coeff = 0.0 * kilocalorie_per_mole
        direction = 1.0

        #soft-core parameters are fixed (the same in all states)
        umsc = float(self.keywords.get('UMAX')) * kilocalorie_per_mole
        ubcore = self.keywords.get('UBCORE')
        if ubcore:
            ubcore = float(ubcore) * kilocalorie_per_mole
        else:
            ubcore = 0.0 * kilocalorie_per_mole
        acore = float(self.keywords.get('ACORE'))

        #create ATM Force
        self.atm_utils.setNonbondedForceGroup(self.nonbondedforcegroup)
        atmvariableforcegroups = [self.nonbondedforcegroup]
        self.atmforce = ATMMetaForce(lambda1, lambda2,  alpha * kilojoules_per_mole, u0/kilojoules_per_mole, w0coeff/kilojoules_per_mole, umsc/kilojoules_per_mole, ubcore/kilojoules_per_mole, acore, direction, atmvariableforcegroups)

        for i in range(self.topology.getNumAtoms()):
            self.atmforce.addParticle(i, 0., 0., 0.)
        for i in self.lig1_atoms:
            self.atmforce.setParticleParameters(i, i, self.displ[0], self.displ[1], self.displ[2] )
        for i in self.lig2_atoms:
            self.atmforce.setParticleParameters(i, i, -self.displ[0], -self.displ[1], -self.displ[2] )

        self.atmforce.setForceGroup(self.atmforcegroup)
        self.system.addForce(self.atmforce)

        #these are the global parameters specified in the cntl files that need to be reset after reading the first configuration
        self.cparams["ATMUmax"] = umsc/kilojoules_per_mole
        self.cparams["ATMUbcore"] = ubcore/kilojoules_per_mole
        self.cparams["ATMAcore"] = acore

    def create_system(self):

        self.load_amber_system()
        self.atm_utils = ATMMetaForceUtils(self.system)
        self.set_ligand_atoms()
        self.set_displacement()
        self.set_vsite_restraints()
        #set orientation restraints
        self.set_orientation_restraints()
        #set reference atoms for alignment force
        self.set_alignmentForce()
        #indexes of the atoms whose position is restrained near the initial positions
        #by a flat-bottom harmonic potential.
        self.set_positional_restraints()
        #temperature is part of the state and is maybe overriden in set_state()
        temperature = 300 * kelvin

        self.set_torsion_metaDbias(temperature)
        
        self.set_atmforce()

        #add barostat
        pressure=1*bar
        self.set_barostat(temperature,pressure,900000000)
        #hack to store ASyncRE quantities in the openmm State
        sforce = mm.CustomBondForce("1")
        for name in self.parameter:
            sforce.addGlobalParameter(self.parameter[name], 0)
        self.system.addForce(sforce)

        self.set_integrator(temperature, self.frictionCoeff, self.MDstepsize)
