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

from async_re import async_re
from bedam_async_re import bedam_async_re_job

from local_openmm_transport import OpenCLContext
from ommreplica import OMMReplica

# OpenMM context overrides for alchemical SDM
class OpenCLContextSDM(OpenCLContext):
    def set_state_values(self, par):
        self.temperature = float(par[0])
        self.lmbd = float(par[1])
        self.lmbd1 = float(par[2])
        self.lmbd2 = float(par[3])
        self.alpha = float(par[4])/kilocalorie_per_mole
        self.u0 = float(par[5]) * kilocalorie_per_mole
        self.w0 = float(par[6]) * kilocalorie_per_mole
        self._inq.put(self.temperature)
        self._inq.put(self.lmbd)
        self._inq.put(self.lmbd1)
        self._inq.put(self.lmbd2)
        self._inq.put(self.alpha)
        self._inq.put(self.u0)
        self._inq.put(self.w0)

    def get_energy_values(self):
        pot_energy = self._outq.get()
        bind_energy = self._outq.get()
        return (pot_energy, bind_energy)

    def _worker_setstate_fromqueue(self):
        temperature = self._inq.get()
        lmbd =  self._inq.get()
        lmbd1 = self._inq.get()
        lmbd2 = self._inq.get()
        alpha = self._inq.get()
        u0 = self._inq.get()
        w0 = self._inq.get()
        self.integrator.setTemperature(temperature)
        self.plugin = self.keywords.get('ATM_PLUGIN')
        if self.plugin == 'ATM-METAFORCE':
            self.atmforce.setLambda1(lmbd1)
            self.atmforce.setLambda2(lmbd2)
            self.atmforce.setAlpha(alpha*kilojoule_per_mole)
            self.atmforce.setU0(u0/ kilojoule_per_mole)
            self.atmforce.setW0(w0 / kilojoule_per_mole)
            self.atmforce.updateParametersInContext(self.context)
        else:
            self.integrator.setLambda(lmbd)
            self.integrator.setLambda1(lmbd1)
            self.integrator.setLambda2(lmbd2)
            self.integrator.setAlpha(alpha*kilojoule_per_mole)
            self.integrator.setU0(u0/ kilojoule_per_mole)
            self.integrator.setW0coeff(w0 / kilojoule_per_mole)
        self.par = [temperature, lmbd, lmbd1, lmbd2, alpha, u0, w0]


    def _worker_getenergy(self):

        if self.plugin == 'ATM-METAFORCE':
            state = self.simulation.context.getState(getEnergy = True, groups = {1,3})
            pot_energy = (state.getPotentialEnergy()).value_in_unit(kilocalorie_per_mole)
            bind_energy = (self.atmforce.getPerturbationEnergy(self.simulation.context)).value_in_unit(kilocalorie_per_mole)
        else:
            bind_energy = (self.integrator.getBindE()*kilojoule_per_mole).value_in_unit(kilocalorie_per_mole)
            pot_energy = (self.integrator.getPotEnergy()*kilojoule_per_mole).value_in_unit(kilocalorie_per_mole)
        self._outq.put(pot_energy)
        self._outq.put(bind_energy)
        self.pot = (pot_energy, bind_energy)

    def  _openmm_worker_body(self):

        self.plugin = self.keywords.get('ATM_PLUGIN')
        if self.plugin == 'ATM-METAFORCE':
            from desmonddmsfile75 import DesmondDMSFile
            from atmmetaforce import ATMMetaForceUtils, ATMMetaForce
        else:
            from simtk.openmm.app.desmonddmsfile import DesmondDMSFile
            from SDMplugin import SDMUtils, LangevinIntegratorSDM
        
        file_input  = '%s_0.dms' % self.basename

        self.dms = DesmondDMSFile(file_input)
        self.topology = self.dms.topology

        self.system = self.dms.createSystem(nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer)

        lig_atoms = self.keywords.get('LIGAND_ATOMS')   #indexes of ligand atoms
        if lig_atoms:
            lig_atoms = [int(i) for i in lig_atoms]
        else:
            msg = "Error: LIGAND_ATOMS is required"
            self._exit(msg)


        if self.plugin == 'ATM-METAFORCE':
            sdm_utils = ATMMetaForceUtils(self.system)
        else:
            sdm_utils = SDMUtils(self.system)
        
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
            cmkf = float(self.keywords.get('CM_KF'))
            kf = cmkf * kilocalorie_per_mole/angstrom**2 #force constant for Vsite CM-CM restraint
            cmtol = float(self.keywords.get('CM_TOL'))
            r0 = cmtol * angstrom #radius of Vsite sphere

            #these can be 'None" if not using orientational restraints
            lig_ref_atoms = None # the 3 atoms of the ligand that define the coordinate system of the ligand
            rcpt_ref_atoms = None # the 3 atoms of the receptor that define the coordinate system of the receptor
            angle_center = None * degrees
            kfangle = None * kilocalorie_per_mole/degrees**2
            angletol = None * degrees
            dihedral1center = None * degrees
            kfdihedral1 = None * kilocalorie_per_mole/degrees**2
            dihedral1tol = None * degrees
            dihedral2center = None * degrees
            kfdihedral2 = None * kilocalorie_per_mole/degrees**2
            dihedral2tol = None * degrees
            
            ligoffset = self.keywords.get('LIGOFFSET')
            if ligoffset:
                ligoffset = [float(offset) for offset in ligoffset.split(',')]*angstrom
            sdm_utils.addRestraintForce(lig_cm_particles = lig_atom_restr,
                                        rcpt_cm_particles = rcpt_atom_restr,
                                        kfcm = kf,
                                        tolcm = r0,
                                        lig_ref_particles = lig_ref_atoms,
                                        rcpt_ref_particles = rcpt_ref_atoms,
                                        angle_center = angle_center,
                                        kfangle = kfangle,
                                        angletol = angletol,
                                        dihedral1center = dihedral1center,
                                        kfdihedral1 = kfdihedral1,
                                        dihedral1tol = dihedral1tol,
                                        dihedral2center = dihedral2center,
                                        kfdihedral2 = kfdihedral2,
                                        dihedral2tol = dihedral2tol,
                                        offset = ligoffset)

        # the integrator object is context-specific
        #temperature = int(self.keywords.get('TEMPERATURES')) * kelvin
        temperature = 300 * kelvin #will be overriden in set_state()
        frictionCoeff = float(self.keywords.get('FRICTION_COEFF')) / picosecond
        MDstepsize = float(self.keywords.get('TIME_STEP')) * picosecond
        umsc = float(self.keywords.get('UMAX')) * kilocalorie_per_mole
        ubcore = self.keywords.get('UBCORE')
        if ubcore:
            ubcore = float(ubcore) * kilocalorie_per_mole
        else:
            ubcore = 0.0 * kilocalorie_per_mole
        acore = float(self.keywords.get('ACORE'))

        if not (self.keywords.get('DISPLACEMENT') is None):
            self.displ = [float(displ) for displ in self.keywords.get('DISPLACEMENT').split(',')]*angstrom
        else:
            msg = "Error: DISPLACEMENT is required"
            self._exit(msg)
        
        if self.plugin == 'ATM-METAFORCE':
            #create ATM Force
            self.atmforce = ATMMetaForce()
        
            for at in self.dms.topology.atoms():
                self.atmforce.addParticle(int(at.id)-1, 0., 0., 0.)
            
            for i in lig_atoms:
                self.atmforce.setParticleParameters(i, i, self.displ[0], self.displ[1], self.displ[2] )

            self.atmforce.setUmax(umsc/kilojoules_per_mole);
            self.atmforce.setUbcore(ubcore/kilojoules_per_mole);
            self.atmforce.setAcore(acore);

            self.atmforce.setForceGroup(3)

            self.system.addForce(self.atmforce)

            self.integrator = LangevinIntegrator(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond )
            self.integrator.setIntegrationForceGroups({1,3})
        else:
            self.integrator = LangevinIntegratorSDM(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond, self.topology.getNumAtoms())
            self.integrator.setBiasMethod(sdm_utils.ILogisticMethod)
            self.integrator.setSoftCoreMethod(sdm_utils.RationalSoftCoreMethod)
            self.integrator.setUmax(umsc / kilojoule_per_mole)
            self.integrator.setAcore(acore)
            self.integrator.setUbcore(ubcore/kilojoule_per_mole)
            for i in lig_atoms:
                self.integrator.setDisplacement(i, self.displ[0]/nanometer, self.displ[1]/nanometer, self.displ[2]/nanometer)
                
class SDMReplica(OMMReplica):
    def __init__(self, replica_id, basename, keywords):
        self.keywords = keywords
        OMMReplica.__init__(self,replica_id, basename)


    #overrides to open dms file for SDM-RE
    def open_dms(self):
        self.plugin = self.keywords.get('ATM_PLUGIN')
        if self.plugin == 'ATM-METAFORCE':
            from desmonddmsfile75 import DesmondDMSFile
        else:
            from simtk.openmm.app.desmonddmsfile import DesmondDMSFile
        
        file_input  = '%s_0.dms' % self.basename

        if not os.path.isdir('r%d' % self._id):
            os.mkdir('r%d' % self._id)

        file_output  = 'r%d/%s_ckp.dms' % (self._id,self.basename)
        if os.path.isfile(file_output):
            print("Reading from %s" % file_output)
        if not os.path.isfile(file_output):
            shutil.copyfile(file_input, file_output)

        self.dms = DesmondDMSFile(file_output) 

        self.sql_conn = self.dms._conn[0]

        # check for sdm_data table in lig dms
        tables = self.dms._tables[0]
        conn = self.sql_conn
        if 'sdm_data' in tables:
            print("Reading sdm_data")
            # read sdm_data table
            q = """SELECT binde,epot,temperature,lambda,lambda1,lambda2,alpha,u0,w0,cycle,stateid,mdsteps FROM sdm_data WHERE id = 1"""
            ans = conn.execute(q)
            for (binde,epot,temperature,lmbd,lambda1,lambda2,alpha,u0,w0,cycle,stateid,mdsteps) in conn.execute(q):
                self.pot = [epot, binde]
                self.par = [temperature, lmbd, lambda1, lambda2, alpha, u0, w0]
                self.cycle = cycle
                self.stateid = stateid
                self.mdsteps = mdsteps
                print("Cycles = %d" % self.cycle)
        else:
            #create sdm_data table with dummy values
            conn.execute("CREATE TABLE IF NOT EXISTS sdm_data (id INTEGER PRIMARY KEY, binde REAL, epot REAL, temperature REAL, lambda REAL, lambda1 REAL, lambda2 REAL, alpha REAL, u0 REAL, w0 REAL, cycle INTEGER, stateid INTEGER, mdsteps INTEGER )")
            conn.execute("INSERT INTO sdm_data (binde,epot,temperature,lambda,lambda1,lambda2,alpha,u0,w0,cycle,stateid,mdsteps) VALUES (0,0,0,0,0,0,0,0,0,0,0,0)")
            conn.commit()
            self.dms._tables[0] = self.dms._readSchemas(conn)

    def save_dms(self):
        self.plugin = self.keywords.get('ATM_PLUGIN')
        if self.plugin == 'ATM-METAFORCE':
            from desmonddmsfile75 import DesmondDMSFile
        else:
            from simtk.openmm.app.desmonddmsfile import DesmondDMSFile
        
        if self.is_state_assigned and self.is_energy_assigned:
            conn = self.sql_conn

            pot_energy =  float(self.pot[0])
            bind_energy = float(self.pot[1])

            temperature = float(self.par[0])
            lmbd =        float(self.par[1])
            lambda1 =     float(self.par[2])
            lambda2 =     float(self.par[3])
            alpha =       float(self.par[4])
            u0 =          float(self.par[5])
            w0coeff =     float(self.par[6])

            conn.execute("UPDATE sdm_data SET binde = %f, epot = %f, temperature = %f, lambda = %f, lambda1 = %f, lambda2 = %f, alpha = %f, u0 = %f, w0 = %f, cycle = %d, stateid = %d, mdsteps = %d WHERE id = 1" % (bind_energy, pot_energy, temperature, lmbd, lambda1, lambda2, alpha, u0, w0coeff, self.cycle, self.stateid, self.mdsteps))
            conn.commit()
            
            self.dms.setPositions(self.positions)                    
            self.dms.setVelocities(self.velocities)

    def set_statepot_from_outputfile(self, replica, cycle):
        outfile = "r%d/%s_%d.out" % (replica, self.basename, cycle)
        data = self._getOpenMMData(outfile)
        # example of format  (ilogistic potential):
        # <temperature, lambda, lambda1, lambda2, alpha, u0, w0,  potential energy, binding energy>
        #   0         1       2       3     4    5       6                7       8
        #
        nr = len(data)
        temperature = float(data[nr-1][0])
        lmbd =    float(data[nr-1][1])
        lambda1 = float(data[nr-1][2])
        lambda2 = float(data[nr-1][3])
        alpha =   float(data[nr-1][4])
        u0    =   float(data[nr-1][5])
        w0    =   float(data[nr-1][6])
        parameters = [temperature,lmbd, lambda1, lambda2, alpha, u0, w0]
        pot_energy = data[nr-1][7]
        binding_energy = data[nr-1][8]
        self.set_state(self.stateid, parameters)
        self.set_energy([pot_energy, binding_energy])

    def set_posvel_from_file(self, replica, cycle):
        self.plugin = self.keywords.get('ATM_PLUGIN')
        if self.plugin == 'ATM-METAFORCE':
            from desmonddmsfile75 import DesmondDMSFile
        else:
            from simtk.openmm.app.desmonddmsfile import DesmondDMSFile
        
        tfile = "r%d/%s_%d.dms" % (replica, self.basename, cycle)
        dms = DesmondDMSFile(tfile)        
        self.positions = copy.deepcopy(dms.positions)
        self.velocities = copy.deepcopy(dms.velocities)
        dms.close()

    def write_posvel_to_file(self, replica, cycle):
        self.plugin = self.keywords.get('ATM_PLUGIN')
        if self.plugin == 'ATM-METAFORCE':
            from desmonddmsfile75 import DesmondDMSFile
        else:
            from simtk.openmm.app.desmonddmsfile import DesmondDMSFile
        
        tfile = "r%d/%s_%d.dms" % (replica, self.basename, cycle)
        dms = DesmondDMSFile(tfile)
        dms.setPositions(self.positions)                    
        dms.setVelocities(self.velocities)
        dms.close()

    def save_out(self):
        pot_energy = self.pot[0]
        bind_energy = self.pot[1]
        temperature = self.par[0]
        lmbd = self.par[1]
        lmbd1 = self.par[2]
        lmbd2 = self.par[3]
        alpha = self.par[4]
        u0 = self.par[5]
        w0 = self.par[6]
        if self.outfile:
            self.outfile.write("%f %f %f %f %f %f %f %f %f\n" % (temperature, lmbd, lmbd1, lmbd2, alpha, u0, w0, pot_energy, bind_energy))
            self.outfile.flush()
          

class bedamtempt_async_re_job(bedam_async_re_job):
    def _setLogger(self):
        self.logger = logging.getLogger("async_re.bedamtempt_async_re")

    def _checkInput(self):
        async_re._checkInput(self)

        #flag to do calculation with ligand displacement in PBC box
        self.pbcdispl = True

        #make sure BEDAM + TEMPERATURE is wanted
        if self.keywords.get('RE_TYPE') != 'BEDAMTEMPT':
            self._exit("RE_TYPE is not BEDAMTEMPT")
        #BEDAM runs with openmm //edit on 10.15
        if self.keywords.get('ENGINE') != 'OPENMM':
            self._exit("ENGINE is not OPENMM")
        #input files
        self.extfiles = self.keywords.get('ENGINE_INPUT_EXTFILES')
        if not (self.extfiles is None):
            if self.extfiles != '':
                self.extfiles = self.extfiles.split(',')
        #list of lambdas
        if self.keywords.get('LAMBDAS') is None:
            self._exit("LAMBDAS needs to be specified")
        self.lambdas = self.keywords.get('LAMBDAS').split(',')
        #list of temperatures
        if self.keywords.get('TEMPERATURES') is None:
            self._exit("TEMPERATURES needs to be specified")
        self.temperatures = self.keywords.get('TEMPERATURES').split(',')

        #parameters for non-linear potentials
        self.gammas = None
        self.bcoeffs = None
        self.w0coeffs = None

        self.lambda1s = None
        self.lambda2s = None
        self.alphas = None
        self.u0s = None
        self.w0coeffs = None

        #quadbias potential
        if self.keywords.get('GAMMAS') is not None:
            self.gammas = self.keywords.get('GAMMAS').split(',')
            self.bcoeffs = self.keywords.get('BCOEFF').split(',')
            self.w0coeffs = self.keywords.get('W0COEFF').split(',')

        #ilogistic potential
        if self.keywords.get('LAMBDA1') is not None:
            self.lambda1s = self.keywords.get('LAMBDA1').split(',')
            self.lambda2s = self.keywords.get('LAMBDA2').split(',')
            self.alphas = self.keywords.get('ALPHA').split(',')
            self.u0s = self.keywords.get('U0').split(',')
            self.w0coeffs = self.keywords.get('W0COEFF').split(',')

        #build parameters for the lambda/temperatures combined states
        self.nreplicas = self._buildBEDAMStates()
        #executive file's directory
        if self.keywords.get('JOB_TRANSPORT') is 'SSH':
            if self.keywords.get('EXEC_DIRECTORY') is None:
                self._exit("EXEC DIRECTORY needs to be specified")
        #added on 10.19.15
        self.implicitsolvent =  self.keywords.get('IMPLICITSOLVENT')
        self.totalsteps = self.keywords.get('PRODUCTION_STEPS')
        self.jobname = self.keywords.get('ENGINE_INPUT_BASENAME')
        self.stepgap = self.keywords.get('PRNT_FREQUENCY')

    def _buildBEDAMStates(self):
        self.stateparams = []
        if self.gammas is not None:
            for (lambd,gamma,b,w0) in zip(self.lambdas,self.gammas,self.bcoeffs,self.w0coeffs):
                for tempt in self.temperatures:
                    st = {}
                    st['lambda'] = lambd
                    st['gamma'] = gamma
                    st['bcoeff'] = b
                    st['w0coeff'] = w0
                    st['temperature'] = tempt
                    self.stateparams.append(st)
        elif self.lambda1s is not None:
            for (lambd,lambda1,lambda2,alpha,u0,w0) in zip(self.lambdas,self.lambda1s,self.lambda2s,self.alphas,self.u0s,self.w0coeffs):
                for tempt in self.temperatures:
                    st = {}
                    st['lambda'] = lambd
                    st['lambda1'] = lambda1
                    st['lambda2'] = lambda2
                    st['alpha'] = alpha
                    st['u0'] = u0
                    st['w0coeff'] = w0
                    st['temperature'] = tempt
                    self.stateparams.append(st)
        else:
            for lambd in self.lambdas:
                for tempt in self.temperatures:
                    st = {}
                    st['lambda'] = lambd
                    st['temperature'] = tempt
                    self.stateparams.append(st)

        return len(self.stateparams)


    def _buildInpFile(self, replica):
        """
        Builds input file for a BEDAM replica based on template input file
        BASENAME.inp for the specified replica at lambda=lambda[stateid] for the
        specified cycle.
        """

        if self.transport_mechanism == "LOCAL_OPENMM":
            return

        basename = self.basename
        stateid = self.status[replica]['stateid_current']
        cycle = self.status[replica]['cycle_current']

        template = "%s.py" % basename
        inpfile = "r%d/%s_%d.py" % (replica, basename, cycle)
        implicitsolvent = self.implicitsolvent
        totalsteps = self.totalsteps
        jobname = self.jobname
        stepgap = self.stepgap


        lambd = self.stateparams[stateid]['lambda']
        temperature = self.stateparams[stateid]['temperature']
        # read template buffer
        tfile = self._openfile(template, "r")
        tbuffer = tfile.read()
        tfile.close()
        # make modifications
        tbuffer = tbuffer.replace("@n@",str(cycle))
        tbuffer = tbuffer.replace("@nm1@",str(cycle-1))
        tbuffer = tbuffer.replace("@lambda@",lambd)
        tbuffer = tbuffer.replace("@temperature@",temperature)
        if 'gamma' in self.stateparams[stateid].keys():
            gamma = self.stateparams[stateid]['gamma']
            b = self.stateparams[stateid]['bcoeff']
            w0 = self.stateparams[stateid]['w0coeff']
            tbuffer = tbuffer.replace("@gamma@",gamma)
            tbuffer = tbuffer.replace("@bcoeff@",b)
            tbuffer = tbuffer.replace("@w0coeff@",w0)

        if 'lambda1' in self.stateparams[stateid].keys():
            lambda1 = self.stateparams[stateid]['lambda1']
            lambda2 = self.stateparams[stateid]['lambda2']
            alpha = self.stateparams[stateid]['alpha']
            u0 = self.stateparams[stateid]['u0']
            w0 = self.stateparams[stateid]['w0coeff']
            tbuffer = tbuffer.replace("@lambda1@",lambda1)
            tbuffer = tbuffer.replace("@lambda2@",lambda2)
            tbuffer = tbuffer.replace("@alpha@",alpha)
            tbuffer = tbuffer.replace("@u0@",u0)
            tbuffer = tbuffer.replace("@w0coeff@",w0)

        # write out
        ofile = self._openfile(inpfile, "w")
        ofile.write(tbuffer)
        ofile.close()

        # update the history status file
        ofile = self._openfile("r%d/state.history" % replica, "a")
        ofile.write("%d %d %s %s\n" % (cycle, stateid, lambd, temperature))
        ofile.close()

    def _extractLast_lambda_BindingEnergy_PotEnergy(self,repl,cycle):
        """
        Extracts binding energy etc. from replica objects
        works only for iLog potential
        """
        replica = self.openmm_replicas[repl]
        (stateid, par) = replica.get_state()
        pot = replica.get_energy()

        pot_energy =  pot[0]
        bind_energy = pot[1]

        temperature = par[0]
        lmbd =        par[1]
        lambda1 =     par[2]
        lambda2 =     par[3]
        alpha =       par[4]
        u0 =          par[5]
        w0 =          par[6]

        if bind_energy == None:
            msg = "Error retrieving state for replica %d" % repl
            self._exit(msg)
        parameters = [temperature, lmbd, lambda1, lambda2, alpha, u0, w0]
        return (parameters, bind_energy, pot_energy)

    def print_status(self):
        """
        Writes to BASENAME_stat.txt a text version of the status of the RE job

        It's fun to follow the progress in real time by doing:

        watch cat BASENAME_stat.txt
        """
        logfile = "%s_stat.txt" % self.basename
        ofile = self._openfile(logfile,"w")
        log = "Replica  State  Lambda Lambda1 Lambda2 Alpha U0 W0coeff Temperature Status  Cycle \n"
        for k in range(self.nreplicas):
            stateid = self.status[k]['stateid_current']
            log += "%6d   %5d  %s %s %s %s %s %s %s %5s  %5d \n" % (k, stateid, self.stateparams[stateid]['lambda'], self.stateparams[stateid]['lambda1'], self.stateparams[stateid]['lambda2'], self.stateparams[stateid]['alpha'], self.stateparams[stateid]['u0'], self.stateparams[stateid]['w0coeff'], self.stateparams[stateid]['temperature'], self.status[k]['running_status'], self.status[k]['cycle_current'])
        log += "Running = %d\n" % self.running
        log += "Waiting = %d\n" % self.waiting

        ofile.write(log)
        ofile.close()

    def _getPot(self,repl,cycle):
        (parameters, u, epot) = self._extractLast_lambda_BindingEnergy_PotEnergy(repl,cycle)

        temperature = float(parameters[0])
        lmb = float(parameters[1])
        lambda1 = float(parameters[2])
        lambda2 = float(parameters[3])
        alpha = float(parameters[4])
        u0 = float(parameters[5])
        w0 = float(parameters[6])

        uf = float(u)
        ee = 1.0 + math.exp(-alpha*(uf-u0))
        ebias = 0.0
        if alpha > 0:
            ebias = ((lambda2 - lambda1)/alpha) * math.log(ee)
        ebias += lambda2 * uf + w0
        e0 = float(epot) - ebias
        return (e0,uf)

    def _getPar(self,repl):
        sid = self.status[repl]['stateid_current']
        lmb = float(self.stateparams[sid]['lambda'])
        tempt = float(self.stateparams[sid]['temperature'])
        kb = 0.0019872041
        beta = 1./(kb*tempt)
        parameters = [beta,lmb]
        lambda1 = float(self.stateparams[sid]['lambda1'])
        lambda2 = float(self.stateparams[sid]['lambda2'])
        alpha =  float(self.stateparams[sid]['alpha'])
        u0 = float(self.stateparams[sid]['u0'])
        w0 = float(self.stateparams[sid]['w0coeff'])
        parameters.append(lambda1)
        parameters.append(lambda2)
        parameters.append(alpha)
        parameters.append(u0)
        parameters.append(w0)

        return (parameters)

    def _reduced_energy(self,par,pot):
        # par: list of parameters
        # pot: list of potentials
        # This is for the ilogistic potential
        beta = par[0]
        lmb = par[1]
        lambda1 = par[2]
        lambda2 = par[3]
        alpha = par[4]
        u0 = par[5]
        w0 = par[6]

        e0 = pot[0]
        uf = pot[1]

        ee = 1.0 + math.exp(-alpha*(uf-u0))
        ebias = 0.0
        if alpha > 0:
            ebias = ((lambda2 - lambda1)/alpha) * math.log(ee)
        ebias += lambda2 * uf + w0
        return beta*(e0 + ebias)

    def checkpointJob(self):
        #disable ctrl-c
        s = signal.signal(signal.SIGINT, signal.SIG_IGN)
        # update replica objects of waiting replicas
        for repl in [k for k in range(self.nreplicas) if self.status[k]['running_status'] == 'W']:
            stateid = self.status[repl]['stateid_current']
            lambd = self.stateparams[stateid]['lambda']
            temperature = float(self.stateparams[stateid]['temperature'])
            lambda1 = float(self.stateparams[stateid]['lambda1'])
            lambda2 = float(self.stateparams[stateid]['lambda2'])
            alpha = float(self.stateparams[stateid]['alpha'])
            u0 = float(self.stateparams[stateid]['u0'])
            w0 = float(self.stateparams[stateid]['w0coeff'])
            par = [temperature, lambd, lambda1, lambda2, alpha, u0, w0]
            self.openmm_replicas[repl].set_state(stateid, par)
        for replica in self.openmm_replicas:
            replica.save_dms()
        signal.signal(signal.SIGINT, s)

    #override for creating SDM versions of the contexts
    def CreateOpenCLContext(self,basename, platform_id = None, device_id = None):
        return OpenCLContextSDM(basename, platform_id, device_id, self.keywords)

    #override for creating SDM versions of the replicas
    def CreateReplica(self, repl_id, basename):
        return SDMReplica(repl_id, basename,self.keywords)

    #override for launching an SDM replica
    def _launchReplica(self,replica,cycle):
        """
        Launches a SDM OpenMM sub-job
        """
        input_file = "%s_%d.py" % (self.basename, cycle)
        log_file = "%s_%d.log" % (self.basename, cycle)
        err_file = "%s_%d.err" % (self.basename, cycle)

        if self.transport_mechanism == "SSH":
            rstfile_rcpt_p = "%s_rcpt_%d.dms" % (self.basename,cycle-1)
            rstfile_lig_p = "%s_lig_%d.dms" % (self.basename,cycle-1)
            local_working_directory = os.getcwd() + "/r" + str(replica)
            remote_replica_dir = "%s_r%d_c%d" % (self.basename, replica, cycle)
            executable = "./runopenmm"

            #sync positions/velocities from internal replica to input dms file
            self.openmm_replicas[replica].write_posvel_to_file(replica, cycle-1)

            job_info = {
                "replica": replica,
                "cycle": cycle,
                "executable": executable,
                "input_file": input_file,
                "output_file": log_file,
                "error_file": err_file,
                "working_directory": local_working_directory,
                "remote_replica_dir": remote_replica_dir,
                "job_input_files": None,
                "job_output_files": None,
                "exec_directory": None}

            # detect which kind of architecture the node use, then choosing
            # different library files and binary files in different lib and bin
            # folders
            if self.keywords.get('EXEC_DIRECTORY'):
                exec_directory = self.keywords.get('EXEC_DIRECTORY')
            else:
                exec_directory = os.getcwd()

            job_info["exec_directory"]=exec_directory

            job_input_files = []
            job_input_files.append(input_file)
            if rstfile_rcpt_p and rstfile_lig_p:
                job_input_files.append(rstfile_rcpt_p)
                job_input_files.append(rstfile_lig_p)
            for filename in self.extfiles:
                job_input_files.append(filename)


            job_output_files = []
            job_output_files.append(log_file)
            job_output_files.append(err_file)
            output_file = "%s_%d.out" % (self.basename, cycle)

            rcptfile="%s_rcpt_%d.dms" % (self.basename,cycle)
            ligfile="%s_lig_%d.dms" % (self.basename,cycle)
            pdbfile="%s_%d.pdb" % (self.basename,cycle)
            dcdfile="%s_%d.dcd" % (self.basename,cycle)

            job_output_files.append(output_file)

            job_output_files.append(rcptfile)
            job_output_files.append(ligfile)
            job_output_files.append(pdbfile)
            job_output_files.append(dcdfile)

            job_info["job_input_files"] = job_input_files
            job_info["job_output_files"] = job_output_files

        elif self.transport_mechanism == "LOCAL_OPENMM":

            nsteps = int(self.keywords.get('PRODUCTION_STEPS'))
            nprnt = int(self.keywords.get('PRNT_FREQUENCY'))
            ntrj = int(self.keywords.get('TRJ_FREQUENCY'))
            if not nprnt % nsteps == 0:
                self._exit("nprnt must be an integer multiple of nsteps.")
            if not ntrj % nsteps == 0:
                sys._exit("ntrj must be an integer multiple of nsteps.")

            job_info = {
                "replica": replica,
                "cycle": cycle,
                "nsteps": nsteps,
                "nprnt": nprnt,
                "ntrj": ntrj
            }

            if self.keywords.get('HEAT_AND_COOL_RATE') is not None:
                probht = float(self.keywords.get('HEAT_AND_COOL_RATE'))
                if probht > random.random():
                    job_info['nheating'] = int(self.keywords.get('HEATING_STEPS'))
                    job_info['ncooling'] = int(self.keywords.get('COOLING_STEPS'))
                    job_info['hightemp'] = float(self.keywords.get('HIGHTEMPERATURE'))

        else: #local with runopenmm?
            executable = os.getcwd() + "/runopenmm" #edit on 10.19
            working_directory = os.getcwd() + "/r" + str(replica)
            job_info = {"executable": executable,
                        "input_file": input_file,
                        "output_file": log_file,
                        "error_file": err_file,
                        "working_directory": working_directory,
                        "cycle": cycle}
            #delete failed file if present
            failed_file = "r%s/%s_%d.failed" % (str(replica),self.basename,cycle)
            if os.path.exists(failed_file):
                os.remove(failed_file)

        if self.keywords.get('VERBOSE') == "yes":
            msg = "_launchReplica(): Launching %s %s in directory %s cycle %d"
            if self.transport_mechanism is 'SSH':
                self.logger.info(msg, executable, input_file, local_working_directory, cycle)
            elif not (self.transport_mechanism ==  'LOCAL_OPENMM'):
                self.logger.info(msg, executable, input_file, working_directory, cycle)

        status = self.transport.launchJob(replica, job_info)

        return status

    def update_state_of_replica(self, repl):
        replica = self.openmm_replicas[repl]

        #retrieve previous state if set
        (old_stateid, old_par) =  replica.get_state()
        if old_stateid != None:
            old_temperature = old_par[0]

        stateid = self.status[repl]['stateid_current']
        temperature = float(self.stateparams[stateid]['temperature'])
        lambd = float(self.stateparams[stateid]['lambda'])
        lambda1 = float(self.stateparams[stateid]['lambda1'])
        lambda2 = float(self.stateparams[stateid]['lambda2'])
        alpha = float(self.stateparams[stateid]['alpha'])
        u0 = float(self.stateparams[stateid]['u0'])
        w0 = float(self.stateparams[stateid]['w0coeff'])
        par = [temperature, lambd, lambda1, lambda2, alpha, u0, w0]
        replica.set_state(stateid, par)

        #rescale velocities (relevant only if state has changed)
        if old_stateid != None:
            if stateid != old_stateid:
                scale = math.sqrt(float(temperature)/float(old_temperature))
                for i in range(0,len(replica.velocities)):
                    replica.velocities[i] = scale*replica.velocities[i]

    def _hasCompleted(self,replica,cycle):
        """
        Returns true if an OpenMM replica has successfully completed a cycle.
        """
        if self.transport_mechanism == "LOCAL_OPENMM":
            #safeguards are off for local transport
            return True

        output_file = "r%s/%s_%d.out" % (replica,self.basename,cycle)
        failed_file = "r%s/%s_%d.failed" % (replica,self.basename,cycle)
        dmsfile_lig = "r%s/%s_lig_%d.dms" % (replica,self.basename,cycle)
        dmsfile_rcpt = "r%s/%s_rcpt_%d.dms" % (replica,self.basename,cycle)

        if os.path.exists(failed_file):
            return False

        #check existence of dms files
        if not self.transport_mechanism == "LOCAL_OPENMM":
            try:
                if not (os.path.exists(dmsfile_rcpt) and os.path.exists(dmsfile_lig)):
                    self.logger.warning("Cannot find file %s and %s", dmsfile_rcpt, dmsfile_lig)
                    return False
            except:
                self.logger.error("Error accessing file %s and %s", dmsfile_rcpt, dmsfile_lig)
                return False

        #check that we can read data from .out
        try:
            self.openmm_replicas[replica]._getOpenMMData(output_file)
        except:
            self.logger.warning("Unable to read/parse output file for replica %d cycle %d" % (replica, cycle))
            return False

        return True


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

    rx = bedamtempt_async_re_job(commandFile, options=None)

    rx.setupJob()

    rx.scheduleJobs()
