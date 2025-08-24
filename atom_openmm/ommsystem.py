from __future__ import print_function
from __future__ import division
"""
Collects all of the ways that openmm systems are loaded
"""
import os, re, sys, time, shutil, copy, random, signal, copy
import numpy as np
import multiprocessing as mp
#from multiprocessing import Process, Queue, Event
import logging

import openmm as mm
from openmm.app import *
from openmm import *
from openmm.unit import *
from datetime import datetime

from atom_openmm.utils.AtomUtils import AtomUtils, separate14

# OpenMM's MTSLangevinIntegrator does not have a setTemperature method
class ATMMTSLangevinIntegrator(MTSLangevinIntegrator):
    def setTemperature(self, temperature):
        self.setGlobalVariableByName('kT', MOLAR_GAS_CONSTANT_R*temperature)

class OMMSystem(object):
    def __init__(self, basename, keywords, pdbtopfile, systemfile, logger):
        self.system = None
        self.topology = None
        self.positions = None
        self.boxvectors = None
        self.integrator = None
        self.barostat = None
        self.keywords = keywords
        self.basename = basename
        self.logger = logger
        self.pdbtopfile = pdbtopfile
        self.systemfile = systemfile
        
        #parameters stored in the openmm state
        self.parameter = {}
        self.parameter['stateid'] = 'REStateId'
        self.parameter['cycle'] = 'RECycle'
        self.parameter['mdsteps'] = 'REMDSteps'
        self.parameter['temperature'] = 'RETemperature'
        self.parameter['potential_energy'] = 'REPotEnergy'
        #more ATM property names are in atmmetaforce

        #parameters from the cntl file
        self.cparams = {}

        self.atmforcegroup = None
        self.nonbondedforcegroup = None
        self.metaDforcegroup = None
        
        self.var_regions = False
        self.var_force_group = None
        if keywords.get('VARIABLE_FORCE_GROUP') is not None:
            self.var_force_group = int(keywords.get('VARIABLE_FORCE_GROUP'))

        self.frictionCoeff = float(self.keywords.get('FRICTION_COEFF')) / picosecond
        self.MDstepsize = float(self.keywords.get('TIME_STEP')) * picosecond

        self.doMetaD = False

        self.major_ommversion = None
        self.minor_ommversion = None
        self.patch_ommversion = None
        version = mm.Platform.getOpenMMVersion().split('.')
        self.logger.info("OpenMM version %s" %  mm.Platform.getOpenMMVersion())
        self.major_ommversion = int(version[0])
        self.minor_ommversion = int(version[1])
        try:
           self.patch_ommversion = int(version[2])
        except:
            pass

        self.v82plus = False
        if (self.major_ommversion > 8) or ( (self.major_ommversion == 8) and (self.minor_ommversion > 1)):
            self.v82plus = True

    def _exit(self, message):
        """Print and flush a message to stdout and then exit."""
        self.logger.error(message)
        sys.stdout.flush()
        sys.exit(1)

    def load_system(self):
        """
        load the topology from a pdb file and the system from an xml file
        """
        self.pdb = PDBFile(self.pdbtopfile)
        self.topology = self.pdb.topology
        self.positions = self.pdb.positions
        self.boxvectors = self.topology.getPeriodicBoxVectors()
        
        #hack to initialize AGBNP's deserializer if present
        try:
            from AGBNPplugin import AGBNPForce
            gb = AGBNPForce()
        except ImportError:
            pass

        #HMASS is set in the system file
        with open(self.systemfile) as input:
            self.system = XmlSerializer.deserialize(input.read())

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

    def free_force_group(self):
        varset = set([self.var_force_group]) if self.var_force_group is not None else set()
        freeGroups = set(range(32)) - set(force.getForceGroup() for force in self.system.getForces()) - varset
        if len(freeGroups) == 0:
            self._exit('Cannot find a free force group. '
                       'The maximum number (32) of the force groups is already used.')
        return max(freeGroups)
        
    def set_positional_restraints(self):
        #indexes of the atoms whose position is restrained near the initial positions
        #by a flat-bottom harmonic potential. 
        posrestr_atoms = self.keywords.get('POS_RESTRAINED_ATOMS')
        self.posrestrForce = None
        if posrestr_atoms is not None:
            fc = float(self.keywords.get('POSRE_FORCE_CONSTANT')) * (kilocalorie_per_mole/angstrom**2)
            tol = float(self.keywords.get('POSRE_TOLERANCE')) * angstrom
            self.posrestrForce = self.atm_utils.addPosRestraints(posrestr_atoms, self.positions, fc, tol)

    def set_torsion_metaDbias(self,temperature):
        from atom_openmm.utils.config import parse_config

        if self.keywords.get('METADBIAS_DIR') == None:
            return
        bias_dirs = self.keywords.get('METADBIAS_DIR')
        bias_offsets = self.keywords.get('METADBIAS_IDXOFFSET')

        for mdir,offset in zip(bias_dirs,bias_offsets) :
            for ext in ("cntl", "yaml", "json"):
                cntlfile = "%s/%s.%s" % (mdir, mdir, ext)
                if os.path.exists(cntlfile):
                    break
            else:
                self._exit("Error: No cntl/yaml/json config file found in %s" % mdir)
            keywords  = parse_config(cntlfile)

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
            self.metaDforcegroup = metaD._force.getForceGroup()
        self.doMetaD = True

    def add_forces_to_atmforce(self):
        import re
        nbpattern = re.compile(".*Nonbonded.*")
        gbpattern = re.compile(".*GB.*")
        harmpattern = re.compile(".*Harmonic.*")
        torpattern  = re.compile(".*Torsion.*")
        if self.var_force_group is not None:
            #add all forces in the var_force_group to ATMForce
            force_removed = True
            while force_removed:
                force_removed = False
                for i in range(self.system.getNumForces()):
                    if self.system.getForce(i).getForceGroup() == self.var_force_group:
                        self.logger.info("Adding Force %s to ATMForce" % self.system.getForce(i).getName())
                        self.atmforce.addForce(copy.copy(self.system.getForce(i)))
                        self.system.removeForce(i)
                        force_removed = True
                        break
        elif self.var_regions:
            #add all standard bonded and non-bonded forces to ATMForce
            force_removed = True
            while force_removed:
                force_removed = False
                for i in range(self.system.getNumForces()):
                    if ( nbpattern.match(self.system.getForce(i).getName())   or
                         gbpattern.match(self.system.getForce(i).getName())   or
                         harmpattern.match(self.system.getForce(i).getName()) or
                         torpattern.match(self.system.getForce(i).getName())    ):
                        self.logger.info("Adding Force %s to ATMForce" % self.system.getForce(i).getName())
                        self.atmforce.addForce(copy.copy(self.system.getForce(i)))
                        self.system.removeForce(i)
                        force_removed = True
                        break
        else:
            #add only non-bonded forces, after separating out the 1-4 interactions
            force_removed = True
            while force_removed:
                force_removed = False
                for i in range(self.system.getNumForces()):
                    if nbpattern.match(self.system.getForce(i).getName()):
                        #separate 1-4 interactions from non-bonded force so they get evaluated
                        #with the bonded terms
                        self.logger.info("Adding Force %s to ATMForce" % self.system.getForce(i).getName())
                        force14 = separate14(self.system.getForce(i))
                        self.system.addForce(force14)
                        self.atmforce.addForce(copy.copy(self.system.getForce(i)))
                        force_removed = True
                        self.system.removeForce(i)
                        break
                    elif gbpattern.match(self.system.getForce(i).getName()):
                        self.logger.info("Adding Force %s to ATMForce" % self.system.getForce(i).getName())
                        self.atmforce.addForce(copy.copy(self.system.getForce(i)))
                        force_removed = True
                        self.system.removeForce(i)
                        break

        self.logger.info("System's Forces after adding variable Forces to ATMForce:")
        for i in range(self.system.getNumForces()):
            self.logger.info("   %s" % self.system.getForce(i).getName())

    def set_integrator(self, temperature, frictionCoeff, MDstepsize, defaultMDstepsize = 0.001*picosecond):
        integrator = self.keywords.get("INTEGRATOR", "atmmts")

        # the Drude polarizable force field requires a special integrator
        drudeForce = None
        # check for the Drude Force in the system
        for i in range(self.system.getNumForces()):
            if self.system.getForce(i).getName() == "DrudeForce":
                drudeForce = copy.copy(self.system.getForce(i))
                break
        # check for the Drude Force in ATMForce
        if drudeForce is None:
            if self.atmforce is not None:
                for i in range(self.atmforce.getNumForces()):
                    if self.atmforce.getForce(i).getName() == "DrudeForce":
                        drudeForce = copy.copy(self.atmforce.getForce(i))
                        break
        if drudeForce is not None:
            integrator = "DrudeLangevin"

        self.logger.info(f"Integrator: {integrator}")

        if integrator.lower() == "atmmts":
            #set the multiplicity of the calculation of bonded forces so that they are evaluated at least once every 1 fs (default time-step)
            bonded_frequency = max(1, int(round(MDstepsize/defaultMDstepsize)))
            self.logger.info("Running with a %f fs time-step with bonded forces integrated %d times per time-step" % (MDstepsize/femtosecond, bonded_frequency))
            if self.doMetaD:
                fgroups = [(0,bonded_frequency), (self.metaDforcegroup, bonded_frequency), (self.atmforcegroup,1)]
            else:
                fgroups = [(0,bonded_frequency), (self.atmforcegroup,1)]
            self.integrator = ATMMTSLangevinIntegrator(temperature, frictionCoeff, MDstepsize, fgroups )
            self.integrator.setConstraintTolerance(0.00001)
        elif integrator.lower() == "langevinmiddle":
            self.integrator = LangevinMiddleIntegrator(
                temperature, frictionCoeff, MDstepsize
            )
            self.integrator.setIntegrationForceGroups({0, self.atmforcegroup})

            self.logger.info(f"Temperature: {self.integrator.getTemperature()}")
            self.logger.info(f"Friction: {self.integrator.getFriction()}")
            self.logger.info(f"Step size: {self.integrator.getStepSize()}")
            self.logger.info(
                f"Constraint tolerance: {self.integrator.getConstraintTolerance()}"
            )
        elif integrator.lower() == "drudelangevin":
            self.logger.info("Using Drude Langevin integrator with a %f fs time-step." % (MDstepsize/femtosecond))
            self.drude_temperature = 1.0*kelvin if self.keywords.get('DRUDE_TEMPERATURE') is None else float(self.keywords.get('DRUDE_TEMPERATURE'))*kelvin
            self.drude_frictionCoeff = 20.0/picosecond if self.keywords.get('DRUDE_FRICTIONCOEFF') is None else float(self.keywords.get('DRUDE_FRICTIONCOEFF'))*kelvin
            self.drude_hardwall = 0.02 if self.keywords.get('DRUDE_HARDWALL') is None else float(self.keywords.get('DRUDE_HARDWALL'))
            self.integrator = DrudeLangevinIntegrator(temperature, frictionCoeff, self.drude_temperature, self.drude_frictionCoeff, MDstepsize)
            self.integrator.setMaxDrudeDistance(self.drude_hardwall) # Drude Hardwall
            self.integrator.setDrudeForce(drudeForce)
        else:
            self._exit(f"Invalid integrator: {integrator}")
        
#Temperature RE
class OMMSystemTRE(OMMSystem):
    def set_integrator(self, temperature, frictionCoeff, MDstepsize, defaultMDstepsize = 0.001*picosecond):
        #place non-bonded force in its own group, assume all other bonded forces are in group 0
        nonbonded = [f for f in self.system.getForces() if isinstance(f, NonbondedForce)][0]
        self.nonbondedforcegroup = self.free_force_group()
        nonbonded.setForceGroup(self.nonbondedforcegroup)
        gbpattern = re.compile(".*GB.*")
        for i in range(self.system.getNumForces()):
            if gbpattern.match(str(type(self.system.getForce(i)))):
                self.system.getForce(i).setForceGroup(self.nonbondedforcegroup)
                break
        #set the multiplicity of the calculation of bonded forces so that they are evaluated at least once every 1 fs (default time-step)
        bonded_frequency = max(1, int(round(MDstepsize/defaultMDstepsize)))
        self.logger.info("Running with a %f fs time-step with bonded forces integrated %d times per time-step" % (MDstepsize/femtosecond, bonded_frequency))
        if doMetaD:
            fgroups = [(0,bonded_frequency), (self.metaDforcegroup, bonded_frequency), (self.nonbondedforcegroup,1)]
        else:
            fgroups = [(0,bonded_frequency), (self.nonbondedforcegroup,1)]
        self.integrator = ATMMTSLangevinIntegrator(temperature, frictionCoeff, MDstepsize, fgroups)
        self.integrator.setConstraintTolerance(0.00001)

    def create_system(self):
        self.load_system()
        #the temperature defines the state and will be overriden in set_state()
        temperature = 300 * kelvin
        #set barostat
        self.set_barostat(temperature,1*bar,0)

        #hack to store ASyncRE quantities in the openmm State
        sforce = mm.CustomBondForce("1")
        for name in self.parameter:
            sforce.addGlobalParameter(self.parameter[name], 0)
        self.system.addForce(sforce)

        self.set_integrator(temperature, self.frictionCoeff, self.MDstepsize)

class OMMSystemABFE(OMMSystem):
    def __init__(self, basename, keywords, pdbtopfile, systemfile, logger):
        super().__init__(basename, keywords, pdbtopfile, systemfile, logger)

        self.parameter['perturbation_energy'] = 'REPertEnergy'
        self.parameter['atmintermediate'] = 'REAlchemicalIntermediate'
        self.parameter['bias_energy'] = 'BiasEnergy'
        self.atmforce = None
        self.lig_atoms = None
        self.lig_atoms0 = None
        self.displ = None

    def set_ligand_atoms(self):
        lig_atoms_in = self.keywords.get('LIGAND_ATOMS')   #indexes of ligand atoms
        if lig_atoms_in is not None:
            self.lig_atoms = [int(i) for i in lig_atoms_in]
        else:
            msg = "Error: LIGAND_ATOMS is required"
            self._exit(msg)

        if not (self.keywords.get('LIGAND_ATOMS0') is None):
            self.lig_atoms0 = [int(i) for i in self.keywords.get('LIGAND_ATOMS0')]

    def set_displacement(self):
        #set displacements
        if not (self.keywords.get('DISPLACEMENT') is None):
            self.displ = self.keywords.get('DISPLACEMENT')*angstrom
        else:
            msg = "Error: DISPLACEMENT is required"
            self._exit(msg)

    def set_vsite_restraints(self):
        #CM-CM Vsite restraints
        lig_atom_restr = self.keywords.get('LIGAND_CM_ATOMS')   #indexes of ligand atoms for CM-CM Vsite restraint
        rcpt_atom_restr = self.keywords.get('RCPT_CM_ATOMS')   #indexes of rcpt atoms for CM-CM Vsite restraint

        cmrestraints_present = (rcpt_atom_restr is not None) and (lig_atom_restr is not None)
        self.vsiterestraintForce = None
        if cmrestraints_present:
            cmkf = float(self.keywords.get('CM_KF'))
            kf = cmkf * kilocalorie_per_mole/angstrom**2 #force constant for Vsite CM-CM restraint
            cmtol = float(self.keywords.get('CM_TOL'))
            r0 = cmtol * angstrom #radius of Vsite sphere
            ligoffset = self.keywords.get('LIGOFFSET')
            if ligoffset is not None:
                ligoffset *= angstrom
            self.vsiterestraintForce = self.atm_utils.addVsiteRestraintForceCMCM(lig_cm_particles = lig_atom_restr,
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


    def set_atmforce(self):
        #these define the state and will be overriden in set_state()
        lmbd = 0.0
        lambda1 = lmbd
        lambda2 = lmbd
        alpha = 0.0 / kilocalorie_per_mole
        uh = 0.0 * kilocalorie_per_mole
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

        #perturbation energy offset
        uoffset = 0.0 * kilocalorie_per_mole
        if self.keywords.get('PERTE_OFFSET') is not None:
            uoffset = float(self.keywords.get('PERTE_OFFSET')) * kilocalorie_per_mole

        #create ATM Force
        referencePotExpression = "select(step(Direction), u0, u1) + "
        alchemicalPotExpression = "select(Lambda2-Lambda1 , ((Lambda2-Lambda1)/Alpha)*log(1+exp(-Alpha*(usc-Uh))) + Lambda2*usc + W0, Lambda2*usc + W0);"
        softCoreExpression = "usc = select(Acore, select(step(u-Ubcore), (Umax-Ubcore)*fsc+Ubcore, u), u);" + \
            "fsc = (z^Acore-1)/(z^Acore+1);" + \
            "z = 1 + 2*(y/Acore) + 2*(y/Acore)^2;" + \
            "y = (u-Ubcore)/(Umax-Ubcore);" + \
	    "u = select(step(Direction), 1, -1)*(u1-(u0 + UOffset))"
        self.atmforce = ATMForce(referencePotExpression + alchemicalPotExpression + softCoreExpression)
        self.atmforce.addGlobalParameter("Lambda1", lambda1);
        self.atmforce.addGlobalParameter("Lambda2", lambda2);
        self.atmforce.addGlobalParameter("Alpha", alpha * kilojoules_per_mole);
        self.atmforce.addGlobalParameter("Uh", uh/kilojoules_per_mole);
        self.atmforce.addGlobalParameter("W0", w0coeff/kilojoules_per_mole);
        self.atmforce.addGlobalParameter("Umax", umsc/kilojoules_per_mole);
        self.atmforce.addGlobalParameter("Ubcore", ubcore/kilojoules_per_mole);
        self.atmforce.addGlobalParameter("Acore", acore);
        self.atmforce.addGlobalParameter("Direction", direction);
        self.atmforce.addGlobalParameter("UOffset", uoffset/kilojoules_per_mole);

        #assign a group to ATMForce for multiple time-steps
        self.atmforcegroup = self.free_force_group()
        self.atmforce.setForceGroup(self.atmforcegroup)
        self.system.addForce(self.atmforce)
        self.nonbondedforcegroup = self.free_force_group()

        #adds Forces of the variable group to ATMForce
        self.add_forces_to_atmforce()

        #adds atoms to ATMForce
        for i in range(self.topology.getNumAtoms()):
            self.atmforce.addParticle(Vec3(0., 0., 0.))
        for i in self.lig_atoms:
            self.atmforce.setParticleParameters(i, Vec3(self.displ[0], self.displ[1], self.displ[2])/nanometer )
        if not (self.keywords.get('LIGAND_ATOMS0') is None):
            for i in self.lig_atoms0:
                self.atmforce.setParticleParameters(i, Vec3(self.displ[0], self.displ[1], self.displ[2])/nanometer,
                                                       Vec3(self.displ[0], self.displ[1], self.displ[2])/nanometer)

        #these are the global parameters specified in the cntl files that need to be reset
        #by the worker after reading the first configuration
        self.cparams[self.atmforce.Umax()] = umsc/kilojoules_per_mole
        self.cparams[self.atmforce.Ubcore()] = ubcore/kilojoules_per_mole
        self.cparams[self.atmforce.Acore()] = acore
        self.cparams['UOffset'] = uoffset/kilojoules_per_mole

    def create_system(self):
        self.load_system()

        self.atm_utils = AtomUtils(self.system)

        self.set_ligand_atoms()

        self.set_displacement()

        self.set_vsite_restraints()

        self.set_orientation_restraints()

        self.set_positional_restraints()

        self.set_atmforce()

        #temperature is part of the state and is maybe overriden in set_state()
        temperature = 300 * kelvin

        #add barostat
        pressure=1*bar
        self.set_barostat(temperature,pressure,0)

        #hack to store ASyncRE quantities in the openmm State
        sforce = mm.CustomBondForce("1")
        for name in self.parameter:
            sforce.addGlobalParameter(self.parameter[name], 0)
        self.system.addForce(sforce)

        self.set_integrator(temperature, self.frictionCoeff, self.MDstepsize)


class OMMSystemRBFE(OMMSystem):
    def __init__(self, basename, keywords, pdbtopfile, systemfile, logger):
        super().__init__(basename, keywords, pdbtopfile, systemfile, logger)

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
            self.lig1_atoms = lig1_atoms_in
        else:
            msg = "Error: LIGAND1_ATOMS is required"
            self._exit(msg)
        if lig2_atoms_in is not None:
            self.lig2_atoms = lig2_atoms_in
        else:
            msg = "Error: LIGAND2_ATOMS is required"
            self._exit(msg)

    def set_displacement(self):
        #set displacements and offsets for ligand 1 and ligand 2
        if self.keywords.get('DISPLACEMENT') is not None:
            self.displ = self.keywords.get('DISPLACEMENT')*angstrom
            #offset from VsiteCM
            ligoffset = [0,0,0]*angstrom
            ligoffset_keyword = self.keywords.get('LIGOFFSET')
            if ligoffset_keyword is not None:
                ligoffset = ligoffset_keyword * angstrom
            #the offset for lig1 is the specified LIGOFFSET
            self.lig1offset = ligoffset
            #the offset for lig2 is the displacement + LIGOFFSET
            d2 = [float(d) for d in self.displ/angstrom]
            f2 = [float(f) for f in ligoffset/angstrom ]
            self.lig2offset = [ d2[i]+f2[i] for i in range(3)]*angstrom
        else:
            msg = "DISPLACEMENT is required"
            self._exit(msg)

        #displacements for initial state and final state, if set 
        self.displ0_lig1 = None
        self.displ1_lig1 = None
        self.displ0_lig2 = None
        self.displ1_lig2 = None
        if ((self.keywords.get('DISPL0_LIG1') is not None) and
            (self.keywords.get('DISPL1_LIG1') is not None) and
            (self.keywords.get('DISPL0_LIG2') is not None) and
            (self.keywords.get('DISPL1_LIG2') is not None) ):
            self.displ0_lig1 = self.keywords.get('DISPL0_LIG1')*angstrom
            self.displ1_lig1 = self.keywords.get('DISPL1_LIG1')*angstrom
            self.displ0_lig2 = self.keywords.get('DISPL0_LIG2')*angstrom
            self.displ1_lig2 = self.keywords.get('DISPL1_LIG2')*angstrom

    def set_vsite_restraints(self):
        #ligand 1 Vsite restraint
        lig1_atom_restr = self.keywords.get('LIGAND1_CM_ATOMS')   #indexes of ligand atoms for CM-CM Vsite restraint

        #ligand 2 Vsite restraint
        lig2_atom_restr = self.keywords.get('LIGAND2_CM_ATOMS')   #indexes of ligand atoms for CM-CM Vsite restraint

        #Vsite restraint receptor atoms
        rcpt_atom_restr = self.keywords.get('RCPT_CM_ATOMS')   #indexes of rcpt atoms for CM-CM Vsite restraint
        if rcpt_atom_restr is None:
            rcpt_atom_restr = self.keywords.get('REST_LIGAND_CMREC_ATOMS')

        cmrestraints_present = (rcpt_atom_restr is not None) and (lig1_atom_restr is not None) and (lig2_atom_restr is not None)

        self.vsiterestraintForce1 = None
        self.vsiterestraintForce2 = None
        if cmrestraints_present:
            cmkf = float(self.keywords.get('CM_KF'))
            kf = cmkf * kilocalorie_per_mole/angstrom**2 #force constant for Vsite CM-CM restraint
            cmtol = float(self.keywords.get('CM_TOL'))
            r0 = cmtol * angstrom #radius of Vsite sphere

            #Vsite restraints for ligands 1 and 2
            self.vsiterestraintForce1 = self.atm_utils.addVsiteRestraintForceCMCM(lig_cm_particles = lig1_atom_restr,
                                        rcpt_cm_particles = rcpt_atom_restr,
                                        kfcm = kf,
                                        tolcm = r0,
                                        offset = self.lig1offset)
            self.vsiterestraintForce2 = self.atm_utils.addVsiteRestraintForceCMCM(lig_cm_particles = lig2_atom_restr,
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

        if refatoms1_cntl is None or refatoms2_cntl is None:
            return

        self.refatoms1 = refatoms1_cntl
        lig1_ref_atoms  = [ self.refatoms1[i]+self.lig1_atoms[0] for i in range(3)]
        self.refatoms2 = refatoms2_cntl
        lig2_ref_atoms  = [ self.refatoms2[i]+self.lig2_atoms[0] for i in range(3)]

        #add alignment force
        self.atm_utils.addAlignmentForce(liga_ref_particles = lig1_ref_atoms,
                                    ligb_ref_particles = lig2_ref_atoms,
                                    kfdispl = float(self.keywords.get('ALIGN_KF_SEP'))*kilocalorie_per_mole/angstrom**2,
                                    ktheta = float(self.keywords.get('ALIGN_K_THETA'))*kilocalorie_per_mole,
                                    kpsi = float(self.keywords.get('ALIGN_K_PSI'))*kilocalorie_per_mole,
                                    offset = self.displ)

    def add_atoms_to_atmforce(self):
        #adds atoms to atmforce using standard fixed molecular displacements
        nodispl = Vec3(0., 0., 0.)
        for i in range(self.topology.getNumAtoms()):
            self.atmforce.addParticle( nodispl )

        for i in self.lig1_atoms:
            if (self.displ0_lig1 is not None) and (self.displ1_lig1 is not None):
                self.atmforce.setParticleParameters(i, Vec3(self.displ1_lig1[0], self.displ1_lig1[1], self.displ1_lig1[2])/nanometer,
                                                        Vec3(self.displ0_lig1[0], self.displ0_lig1[1], self.displ0_lig1[2])/nanometer)
            else:
                self.atmforce.setParticleParameters(i,  Vec3(self.displ[0], self.displ[1], self.displ[2])/nanometer, nodispl )

        for i in self.lig2_atoms:
            if (self.displ0_lig2 is not None) and (self.displ1_lig2 is not None):
                self.atmforce.setParticleParameters(i, Vec3(self.displ1_lig2[0], self.displ1_lig2[1], self.displ1_lig2[2])/nanometer,
                                                    Vec3(self.displ0_lig2[0], self.displ0_lig2[1], self.displ0_lig2[2])/nanometer)
            else:
                self.atmforce.setParticleParameters(i, -Vec3(self.displ[0], self.displ[1], self.displ[2])/nanometer, nodispl)

    def add_common_var_atoms_to_atmforce(self):
        #add atoms to atmforce using common/variable regions
        for i in range(self.topology.getNumAtoms()):
            self.atmforce.addParticle( )

        lig1_var_atoms = [int(i) for i in self.keywords.get('LIGAND1_VAR_ATOMS')  ]
        lig2_var_atoms = [int(i) for i in self.keywords.get('LIGAND2_VAR_ATOMS')  ]

        lig1_attach_atom = int(self.keywords.get('LIGAND1_ATTACH_ATOM'))
        lig2_attach_atom = int(self.keywords.get('LIGAND2_ATTACH_ATOM'))

        if self.keywords.get('LIGAND1_COMMON_ATOMS') is None:
            lig1_common_atoms = sorted([i for i in self.lig1_atoms if i not in lig1_var_atoms ])
            lig2_common_atoms = sorted([i for i in self.lig2_atoms if i not in lig2_var_atoms ])
        else:
            lig1_common_atoms = [int(i) for i in self.keywords.get('LIGAND1_COMMON_ATOMS')  ]
            lig2_common_atoms = [int(i) for i in self.keywords.get('LIGAND2_COMMON_ATOMS')  ]

        if not len(lig1_common_atoms) == len(lig2_common_atoms):
            msg = "Error: the number of commong atoms of lig1 (%d) and lig2 (%d) differ" % (len(lig1_common_atoms),len(lig2_common_atoms))
            self._exit(msg)

        try:
            # try official OpenMM>=8.4 ATMForce API
            for i in range(len(lig1_common_atoms)):
                self.atmforce.setParticleTransformation(lig1_common_atoms[i], ParticleOffsetDisplacement(lig2_common_atoms[i], lig1_common_atoms[i]))
            for i in range(len(lig2_common_atoms)):
                self.atmforce.setParticleTransformation(lig2_common_atoms[i], ParticleOffsetDisplacement(lig1_common_atoms[i], lig2_common_atoms[i]))
            for i in range(len(lig1_var_atoms)):
                self.atmforce.setParticleTransformation(lig1_var_atoms[i], ParticleOffsetDisplacement(lig2_attach_atom, lig1_attach_atom))
            for i in range(len(lig2_var_atoms)):
                self.atmforce.setParticleTransformation(lig2_var_atoms[i], ParticleOffsetDisplacement(lig1_attach_atom, lig2_attach_atom))
        except:
            try:
                # try unofficial Gallicchio-Lab atm-coordinate-swap branch OpenMM 8.2 ATMForce API
                for i in range(len(lig1_common_atoms)):
                    self.atmforce.setParticleParameters(lig1_common_atoms[i], lig2_common_atoms[i], lig1_common_atoms[i], -1, -1)
                for i in range(len(lig2_common_atoms)):
                    self.atmforce.setParticleParameters(lig2_common_atoms[i], lig1_common_atoms[i], lig2_common_atoms[i], -1, -1)
                for i in range(len(lig1_var_atoms)):
                    self.atmforce.setParticleParameters(lig1_var_atoms[i], lig2_attach_atom, lig1_attach_atom, -1, -1)
                for i in range(len(lig2_var_atoms)):
                    self.atmforce.setParticleParameters(lig2_var_atoms[i], lig1_attach_atom, lig2_attach_atom, -1, -1)
            except:
                self._exit("Variable displacements are not supported by the OpenMM backend")

    def set_atmforce(self):
        #these define the state and will be overriden in set_state()
        temperature = 300 * kelvin
        lmbd = 0.0
        lambda1 = lmbd
        lambda2 = lmbd
        alpha = 0.0 / kilocalorie_per_mole
        uh = 0.0 * kilocalorie_per_mole
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
        
        #perturbation energy offset
        uoffset = 0.0 * kilocalorie_per_mole
        if self.keywords.get('PERTE_OFFSET') is not None:
            uoffset = float(self.keywords.get('PERTE_OFFSET')) * kilocalorie_per_mole

        #create ATM Force
        referencePotExpression = "select(step(Direction), u0, u1) + "
        alchemicalPotExpression = "select(Lambda2-Lambda1 , ((Lambda2-Lambda1)/Alpha)*log(1+exp(-Alpha*(usc-Uh))) + Lambda2*usc + W0, Lambda2*usc + W0);"
        softCoreExpression = "usc = select(Acore, select(step(u-Ubcore), (Umax-Ubcore)*fsc+Ubcore, u), u);" + \
            "fsc = (z^Acore-1)/(z^Acore+1);" + \
            "z = 1 + 2*(y/Acore) + 2*(y/Acore)^2;" + \
            "y = (u-Ubcore)/(Umax-Ubcore);" + \
	    "u = select(step(Direction), 1, -1)*(u1-(u0 + UOffset))"
        self.atmforce = ATMForce(referencePotExpression + alchemicalPotExpression + softCoreExpression)
        self.atmforce.addGlobalParameter("Lambda1", lambda1);
        self.atmforce.addGlobalParameter("Lambda2", lambda2);
        self.atmforce.addGlobalParameter("Alpha", alpha * kilojoules_per_mole);
        self.atmforce.addGlobalParameter("Uh", uh/kilojoules_per_mole);
        self.atmforce.addGlobalParameter("W0", w0coeff/kilojoules_per_mole);
        self.atmforce.addGlobalParameter("Umax", umsc/kilojoules_per_mole);
        self.atmforce.addGlobalParameter("Ubcore", ubcore/kilojoules_per_mole);
        self.atmforce.addGlobalParameter("Acore", acore);
        self.atmforce.addGlobalParameter("Direction", direction);
        self.atmforce.addGlobalParameter("UOffset", uoffset/kilojoules_per_mole);
        
        #assign a group to ATMForce for multiple time-steps
        self.atmforcegroup = self.free_force_group()
        self.atmforce.setForceGroup(self.atmforcegroup)
        self.system.addForce(self.atmforce)
        self.nonbondedforcegroup = self.free_force_group()

        #displacements based on position of attachment atoms
        self.pos_displacement = False
        if self.keywords.get('LIGAND1_ATTACH_ATOM') is not None:
            self.pos_displacement = True

        #common and variable regions protocol
        self.var_regions = False
        if self.keywords.get('LIGAND1_VAR_ATOMS') is not None:
            self.var_regions = True

        #adds Forces of the given group to ATMForce
        self.add_forces_to_atmforce()

        #adds atoms to ATMForce
        if self.pos_displacement:
            #use common/variable regions
            self.add_common_var_atoms_to_atmforce()
        else:
            #standard displacements
            self.add_atoms_to_atmforce()

        #these are the global parameters specified in the cntl files that need to be reset
        #by the worker after reading the first configuration
        self.cparams[self.atmforce.Umax()] = umsc/kilojoules_per_mole
        self.cparams[self.atmforce.Ubcore()] = ubcore/kilojoules_per_mole
        self.cparams[self.atmforce.Acore()] = acore
        self.cparams['UOffset'] = uoffset/kilojoules_per_mole

    def create_system(self):

        self.load_system()
        self.atm_utils = AtomUtils(self.system)
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
        self.set_barostat(temperature,pressure,0)
        #hack to store ASyncRE quantities in the openmm State
        sforce = mm.CustomBondForce("1")
        for name in self.parameter:
            sforce.addGlobalParameter(self.parameter[name], 0)
        self.system.addForce(sforce)

        self.set_integrator(temperature, self.frictionCoeff, self.MDstepsize)
