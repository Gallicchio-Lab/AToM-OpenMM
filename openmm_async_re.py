import os
import re
import random
import math
import logging
from async_re import async_re

from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from simtk.openmm.app.desmonddmsfile import *
from SDMplugin import *
from local_openmm_transport import OpenCLContext

class openmm_job(async_re):
    def __init__(self, command_file, options):
        async_re.__init__(self, command_file, options)
        
        if self.transport_mechanism == "LOCAL_OPENMM":
            from local_openmm_transport import SDMReplica
            
            # creates openmm context objects
            self.openmm_contexts = []
            pattern = re.compile('(\d+):(\d+)')
            for node in self.compute_nodes:
                slot_id = node["slot_number"]
                matches = pattern.search(slot_id)
                platform_id = int(matches.group(1))
                device_id = int(matches.group(2))
                self.openmm_contexts.append(self.CreateOpenCLContext(self.basename, platform_id, device_id))

            #creates openmm replica objects
            self.openmm_replicas = []
            for i in range(self.nreplicas):
                self.openmm_replicas.append(SDMReplica(i, self.basename))


    def CreateOpenCLContext(self,basename, platform_id = None, device_id = None):
        rcptfile_input  = '%s_rcpt_0.dms' % basename 
        ligfile_input   = '%s_lig_0.dms'  % basename
        
        # the dms object is replica-specific
        dms = DesmondDMSFile([ligfile_input, rcptfile_input]) 
        topology = dms.topology
    
        # the system object is context-specific
        #    system = dms.createSystem(nonbondedMethod=NoCutoff, OPLS = True, implicitSolvent= 'AGBNP')

        #implicit solvent literal not working
        if self.keywords.get('IMPLICITSOLVENT') is not None:
            implicitsolvent = self.keywords.get('IMPLICITSOLVENT')
            #print type(implicitsolvent)
        
        system = dms.createSystem(nonbondedMethod=NoCutoff, OPLS = True, implicitSolvent= None)

        natoms_ligand = int(self.keywords.get('NATOMS_LIGAND'))
        lig_atoms = range(natoms_ligand)
        # atom indexes here refer to indexes in either lig or rcpt dms file, rather than in the complex 
        #lig_atom_restr = [0, 1, 2, 3, 4, 5]   #indexes of ligand atoms for CM-CM Vsite restraint

        cm_lig_atoms = self.keywords.get('REST_LIGAND_CMLIG_ATOMS')   #indexes of ligand atoms for CM-CM Vsite restraint
        #convert the string of lig atoms to integer
        lig_atom_restr = [int(i) for i in cm_lig_atoms]
        #rcpt_atom_restr = [121, 210, 281, 325, 406, 527, 640, 650, 795, 976, 1276]   #indexes of rcpt atoms for CM-CM Vsite restraint
        
        cm_rcpt_atoms = self.keywords.get('REST_LIGAND_CMREC_ATOMS')   #indexes of rcpt atoms for CM-CM Vsite restraint
        #convert the string of receptor rcpt atoms to integer
        rcpt_atom_restr = [int(i) for i in cm_rcpt_atoms]

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

        #transform indexes of receptor atoms
        for i in range(len(rcpt_atom_restr)):
            rcpt_atom_restr[i] += natoms_ligand
            if rcpt_ref_atoms:
                for i in range(len(rcpt_ref_atoms)):
                    rcpt_ref_atoms[i] += natoms_ligand

        sdm_utils = SDMUtils(system, lig_atoms)
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
                                dihedral2tol = dihedral2tol)
    
        # the integrator object is context-specific
        #temperature = int(self.keywords.get('TEMPERATURES')) * kelvin
        temperature = 300 * kelvin
        frictionCoeff = float(self.keywords.get('FRICTION_COEFF')) / picosecond
        MDstepsize = float(self.keywords.get('TIME_STEP')) * picosecond
        umsc = float(self.keywords.get('UMAX')) * kilocalorie_per_mole
        acore = float(self.keywords.get('ACORE'))
        integrator = LangevinIntegratorSDM(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond, lig_atoms)
        integrator.setBiasMethod(sdm_utils.ILogisticMethod)
        integrator.setSoftCoreMethod(sdm_utils.RationalSoftCoreMethod)
        integrator.setUmax(umsc / kilojoule_per_mole)
        integrator.setAcore(acore)

        platform_properties = {}                                 
        platform = Platform.getPlatformByName('OpenCL')
        platformId = platform_id
        deviceId = device_id

        return OpenCLContext(topology, system, integrator, platform, platformId, deviceId, platform_properties)

    
                
    def _setLogger(self):
        self.logger = logging.getLogger("async_re.openmm_async_re")

    def _launchReplica(self,replica,cycle): 
        """
        Launches a OpenMM sub-job
        """
        input_file = "%s_%d.py" % (self.basename, cycle)
        log_file = "%s_%d.log" % (self.basename, cycle)
        err_file = "%s_%d.err" % (self.basename, cycle)
        
        if self.transport_mechanism == "SSH":
	    rstfile_rcpt_p = "%s_rcpt_%d.dms" % (self.basename,cycle-1)
	    rstfile_lig_p = "%s_lig_%d.dms" % (self.basename,cycle-1)
            local_working_directory = os.getcwd() + "/r" + str(replica)
            remote_replica_dir = "%s_r%d_c%d" % (self.basename, replica, cycle)
            executable = "./runopenmm" #edit 10.20

            job_info = {
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
                job_input_files.append(rstfile_rcpt_p) #edit 10.16
	    	job_input_files.append(rstfile_lig_p) #edit 10.16
            for filename in self.extfiles:
                job_input_files.append(filename)


            job_output_files = []
            job_output_files.append(log_file)
            job_output_files.append(err_file)
            output_file = "%s_%d.out" % (self.basename, cycle)

            if self.keywords.get('RE_TYPE') == 'TEMPT':
                dmsfile = "%s_%d.dms" % (self.basename, cycle)
            elif self.keywords.get('RE_TYPE') == 'BEDAMTEMPT':
                rcptfile="%s_rcpt_%d.dms" % (self.basename,cycle)
                ligfile="%s_lig_%d.dms" % (self.basename,cycle)
                pdbfile="%s_%d.pdb" % (self.basename,cycle)
                dcdfile="%s_%d.dcd" % (self.basename,cycle)
                
            job_output_files.append(output_file)

            if self.keywords.get('RE_TYPE') == 'TEMPT':
                job_output_files.append(dmsfile)
            elif self.keywords.get('RE_TYPE') == 'BEDAMTEMPT':
                job_output_files.append(rcptfile)
                job_output_files.append(ligfile)
                job_output_files.append(pdbfile)
                job_output_files.append(dcdfile)
                
            job_info["job_input_files"] = job_input_files;
            job_info["job_output_files"] = job_output_files;

        elif self.transport_mechanism == "LOCAL_OPENMM":
            from local_openmm_transport import SDMReplica
            
            stateid = self.status[replica]['stateid_current']
            lambd = self.stateparams[stateid]['lambda']
            temperature = self.stateparams[stateid]['temperature']
            lambda1 = self.stateparams[stateid]['lambda1']
            lambda2 = self.stateparams[stateid]['lambda2']
            alpha = self.stateparams[stateid]['alpha']
            u0 = self.stateparams[stateid]['u0']
            w0 = self.stateparams[stateid]['w0coeff']
            self.openmm_replicas[replica].set_state(lambd, lambda1, lambda2, alpha, u0, w0)
            nsteps = int(self.keywords.get('PRODUCTION_STEPS'))
            print("DEBUG: setting replica state for replica %d, steps = %d" % (replica, nsteps))
            job_info = {
                "cycle": cycle,
                "nsteps": nsteps
            }
            
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
            elif not self.transport_mechanism is 'LOCAL_OPENMM':
                self.logger.info(msg, executable, input_file, working_directory, cycle)

        status = self.transport.launchJob(replica, job_info)
        
        return status

    def _getOpenMMData(self,file):
        """
        Reads all of the Openmm simulation data values temperature, energies,
        etc.  at each time step and puts into a big table
        """
        if not os.path.exists(file):
            msg = 'File does not exists: %s' % file
            self._exit(msg)
        
        data = []
        f = self._openfile(file, "r")
        line = f.readline()
        while line:
            datablock = []
            for word in line.split():
                datablock.append(float(word))
            data.append(datablock)
            line = f.readline()
        f.close
        return data

    def _hasCompleted(self,replica,cycle):
        """
        Returns true if an OpenMM replica has successfully completed a cycle.
        """
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
        #        try:
        datai = self._getOpenMMData(output_file)
        nf = len(datai[0])
        nr = len(datai)
        #except:
	#    self.logger.warning("Unable to read/parse file %s", output_file)
        #    return False

        return True

    def _computeSwapMatrix(self, replicas, states):
        """
        Compute matrix of dimension-less energies: each column is a replica
        and each row is a state so U[i][j] is the energy of replica j in state
        i.

        Note that the matrix is sized to include all of the replicas and states
        but the energies of replicas not in waiting state, or those of waiting
        replicas for states not belonging to waiting replicas list are
        undefined.
        """
        # U will be sparse matrix, but is convenient bc the indices of the
        # rows and columns will always be the same.
        U = [[ 0. for j in range(self.nreplicas)]
             for i in range(self.nreplicas)]

        n = len(replicas)

        #collect replica parameters and potentials
        par = []
        pot = []
        for k in replicas:
            v = self._getPot(k,self.status[k]['cycle_current'])
            l = self._getPar(k)
            par.append(l)
            pot.append(v)
        if self.verbose:
            self.logger.info("Swap matrix info:")
            self.logger.info("%s", ' '.join(map(str, pot)))
            self.logger.info("%s", ' '.join(map(str, par)))

        for i in range(n):
            repl_i = replicas[i]
            for j in range(n):
                sid_j = states[j]
                # energy of replica i in state j
                U[sid_j][repl_i] = self._reduced_energy(par[j],pot[i])
        return U
