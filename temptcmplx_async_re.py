from __future__ import print_function
from __future__ import division
import os
import sys
import re
import time
import math
import random
import logging
import shutil
from async_re import async_re
from tempt_async_re import tempt_async_re_job
from local_openmm_transport import OpenCLContext
from ommreplica import OMMReplica

from simtk import openmm as mm
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from simtk.openmm.app.desmonddmsfile import *
from datetime import datetime
from SDMplugin import *


# OpenMM context overrides for T-RE for a complex
class OpenCLContextCmplx(OpenCLContext):
    def _openmm_worker_body(self):
        input_dms_file_lig = '%s_lig_0.dms' % self.basename
        input_dms_file_rcpt  = '%s_rcpt_0.dms' % self.basename
        self.dms = DesmondDMSFile([input_dms_file_lig, input_dms_file_rcpt])
        self.topology = self.dms.topology
        implicitsolvent = str(self.keywords.get('IMPLICITSOLVENT'))
        if implicitsolvent is None:
            self.system = self.dms.createSystem(nonbondedMethod=NoCutoff, OPLS = True, implicitSolvent = None)
        elif implicitsolvent == 'AGBNP':
            self.system = self.dms.createSystem(nonbondedMethod=NoCutoff, OPLS = True, implicitSolvent = 'AGBNP')
        else:
            print('Unknown implicit solvent %s' % implicitsolvent)
            sys.exit(1)

        temperature = 300 * kelvin #will be overridden by set_state()
        frictionCoeff = float(self.keywords.get('FRICTION_COEFF')) / picosecond
        MDstepsize = float(self.keywords.get('TIME_STEP')) * picosecond
        self.integrator = LangevinIntegrator(temperature/kelvin, frictionCoeff/(1/picosecond), MDstepsize/ picosecond)

class CmplxReplica(OMMReplica):
    #overrides to open dms file for a complex
    def open_dms(self):
        rcptfile_input  = '%s_rcpt_0.dms' % self.basename
        ligfile_input   = '%s_lig_0.dms'  % self.basename

        if not os.path.isdir('r%d' % self._id):
            os.mkdir('r%d' % self._id)

        ligfile_output  = 'r%d/%s_lig_ckp.dms' % (self._id,self.basename)
        if not os.path.isfile(ligfile_output):
            shutil.copyfile(ligfile_input, ligfile_output)

        rcptfile_output = 'r%d/%s_rcpt_ckp.dms' % (self._id,self.basename)
        if not os.path.isfile(rcptfile_output):
            shutil.copyfile(rcptfile_input, rcptfile_output)

        self.dms = DesmondDMSFile([ligfile_output, rcptfile_output])

        self.sql_conn_lig = self.dms._conn[0]
        self.sql_conn_rcpt = self.dms._conn[1]

        # check for tre_data table in dms file
        tables = self.dms._tables[0]
        conn = self.sql_conn_lig
        if 'tre_data' in tables:
            # read tre_data table
            q = """SELECT epot,temperature,cycle,stateid,mdsteps FROM tre_data WHERE id = 1"""
            ans = conn.execute(q)
            for (epot,temperature,cycle,stateid,mdsteps) in conn.execute(q):
                self.pot = [epot]
                self.par = [temperature]
                self.cycle = cycle
                self.stateid = stateid
                self.mdsteps = mdsteps

    def save_dms(self):
        if self.is_state_assigned and self.is_energy_assigned:
            conn_lig = self.sql_conn_lig
            conn_rcpt = self.sql_conn_rcpt
            conn = conn_lig
            tables = self.dms._tables[0]
            if 'tre_data' not in tables:
                conn.execute("CREATE TABLE IF NOT EXISTS tre_data (id INTEGER PRIMARY KEY, epot REAL, temperature REAL, cycle INTEGER, stateid INTEGER, mdsteps INTEGER )")
                conn.execute("INSERT INTO tre_data (epot,temperature,cycle,stateid,mdsteps) VALUES (0,0,0,0,0)")
                conn.commit()
                self.dms._tables[0] = self.dms._readSchemas(conn)
            pot_energy =  float(self.pot[0])
            temperature = float(self.par[0])
            conn.execute("UPDATE tre_data SET epot = %f, temperature = %f, cycle = %d, stateid = %d, mdsteps = %d WHERE id = 1" % (pot_energy, temperature, self.cycle, self.stateid, self.mdsteps))
            conn.commit()
            self.dms.setPositions(self.positions)
            self.dms.setVelocities(self.velocities)

    def set_posvel_from_file(self, replica, cycle):
        ligfile = "r%d/%s_lig_%d.dms" % (replica, self.basename, cycle)
        rcptfile = "r%d/%s_rcpt_%d.dms" % (replica, self.basename, cycle)
        dms = DesmondDMSFile([ligfile, rcptfile])
        self.positions = copy.deepcopy(dms.positions)
        self.velocities = copy.deepcopy(dms.velocities)
        dms.close()

    def write_posvel_to_file(self, replica, cycle):
        ligfile = "r%d/%s_lig_%d.dms" % (replica, self.basename, cycle)
        rcptfile = "r%d/%s_rcpt_%d.dms" % (replica, self.basename, cycle)
        dms = DesmondDMSFile([ligfile, rcptfile])
        dms.setPositions(self.positions)
        dms.setVelocities(self.velocities)
        dms.close()

class temptcmplx_async_re_job(tempt_async_re_job):
    def CreateOpenCLContext(self,basename, platform_id = None, device_id = None):
        return OpenCLContextCmplx(basename, platform_id, device_id, self.keywords)

    def CreateReplica(self, repl_id, basename):
        return CmplxReplica(repl_id, basename)

if __name__ == '__main__':

    # Parse arguments:
    usage = "%prog <ConfigFile>"

    if len(sys.argv) != 2:
        print("Please specify ONE input file")
        sys.exit(1)

    commandFile = sys.argv[1]

    print("")
    print("====================================")
    print("Temperature Asynchronous Replica Exchange ")
    print("====================================")
    print("")
    print("Started at: " + str(time.asctime()))
    print("Input file:", commandFile)
    print("")
    sys.stdout.flush()

    rx = temptcmplx_async_re_job(commandFile, options=None)

    rx.setupJob()

    rx.scheduleJobs()
