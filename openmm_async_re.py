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
from local_openmm_transport import OMMReplica

class openmm_job(async_re):
    def __init__(self, command_file, options):
        async_re.__init__(self, command_file, options)
        
        if self.transport_mechanism == "LOCAL_OPENMM":
            
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
                replica = self.CreateReplica(i, self.basename)
                self.openmm_replicas.append(replica)

    def CreateOpenCLContext(self,basename, platform_id = None, device_id = None):
        return OpenCLContext(basename, platform_id, device_id, self.keywords)

    def CreateReplica(self, repl_id, basename):
        return OMMReplica(repl_id, basename)
    
    def _setLogger(self):
        self.logger = logging.getLogger("async_re.openmm_async_re")

    # default for temperature replica exchange
    def _launchReplica(self,replica,cycle):
        """
        Launches a T-RE OpenMM sub-job
        """
        input_file = "%s_%d.py" % (self.basename, cycle)
        log_file = "%s_%d.log" % (self.basename, cycle)
        err_file = "%s_%d.err" % (self.basename, cycle)

        if self.transport_mechanism == "SSH":

            # WORK IN PROGRESS
            
            rstfile_p = "%s_%d.dms" % (self.basename,cycle-1)
            local_working_directory = os.getcwd() + "/r" + str(replica)
            remote_replica_dir = "%s_r%d_c%d" % (self.basename, replica, cycle)
            executable = "./runopenmm"

            job_info = {
                "executable": executable,
                "input_file": input_file,
                "output_file": log_file,
                "error_file": err_file,
                "working_directory": local_working_directory,
                "remote_replica_dir": remote_replica_dir,
                "job_input_files": None,
                "job_output_files": None,
                "exec_directory": None
            }

            if self.keywords.get('EXEC_DIRECTORY'):
                exec_directory = self.keywords.get('EXEC_DIRECTORY')
            else:
                exec_directory = os.getcwd()

            job_info["exec_directory"]=exec_directory

            job_input_files = []
            
            job_input_files.append(input_file)
            job_input_files.append(rstfile_p)
            for filename in self.extfiles:
                job_input_files.append(filename)

            job_output_files = []
            job_output_files.append(log_file)
            job_output_files.append(err_file)
            output_file = "%s_%d.out" % (self.basename, cycle)

            dmsfile = "%s_%d.dms" % (self.basename, cycle)
            pdbfile = "%s_%d.pdb" % (self.basename,cycle)
            dcdfile = "%s_%d.dcd" % (self.basename,cycle)
            
            job_output_files.append(output_file)

            job_output_files.append(dmsfile)
            job_output_files.append(pdbfile)
            job_output_files.append(dcdfile)
            
        elif self.transport_mechanism == "LOCAL_OPENMM":
            
            stateid = self.status[replica]['stateid_current']
            temperature = self.stateparams[stateid]['temperature']
            par = (temperature)
            self.openmm_replicas[replica].set_state(stateid, par)
            nsteps = int(self.keywords.get('PRODUCTION_STEPS'))
            job_info = {
                "cycle": cycle,
                "nsteps": nsteps
            }
            
        else: #local with runopenmm?
            self._exit("Unknown job transport")
            
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
