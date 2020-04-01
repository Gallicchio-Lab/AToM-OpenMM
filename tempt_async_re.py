from __future__ import print_function
from __future__ import division
import os
import sys
import re
import time
import math
import random
import logging
from async_re import async_re
from openmm_async_re import openmm_job
from ommreplica import OMMReplica

class tempt_async_re_job(openmm_job):

    def _setLogger(self):
        self.logger = logging.getLogger('async_re.tempt_async_re')


    def _checkInput(self):
        async_re._checkInput(self)
        #make sure TEMPERATURE is wanted
        if self.keywords.get('RE_TYPE') != 'TEMPT':
            self._exit("RE_TYPE is not TEMPT")
        #we run with OPENMM
        if self.keywords.get('ENGINE') != 'OPENMM':
            self._exit("ENGINE is not OPENMM")
        #input files
        self.extfiles = self.keywords.get('ENGINE_INPUT_EXTFILES')
        if not (self.extfiles is None):
            if self.extfiles != '':
                self.extfiles = self.extfiles.split(',')
        #list of temperatures
        if self.keywords.get('TEMPERATURES') is None:
            self._exit("TEMPERATURES needs to be specified")
        temperatures = self.keywords.get('TEMPERATURES').split(',')
        #build parameters for the lambda/temperatures combined states
        self.nreplicas = self._buildStates(temperatures)
        #executive file's directory
        if self.keywords.get('JOB_TRANSPORT') is 'SSH':
            if self.keywords.get('EXEC_DIRECTORY') is None:
                self._exit("EXEC DIRECTORY needs to be specified")

    def _buildStates(self,temperatures):
        self.stateparams = []
        for tempt in temperatures:
            st = {}
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

        temperature = self.stateparams[stateid]['temperature']
        # read template buffer
        tfile = self._openfile(template, "r")
        tbuffer = tfile.read()
        tfile.close()
        # make modifications
        tbuffer = tbuffer.replace("@n@",str(cycle))
        tbuffer = tbuffer.replace("@nm1@",str(cycle-1))
        tbuffer = tbuffer.replace("@temperature@",temperature)
        tbuffer = tbuffer.replace("@jobname@",basename)
        tbuffer = tbuffer.replace("@replica@",str(replica))
        # write out
        ofile = self._openfile(inpfile, "w")
        ofile.write(tbuffer)
        ofile.close()

        # update the history status file
        ofile = self._openfile("r%d/state.history" % replica, "a")
        ofile.write("%d %d %s\n" % (cycle, stateid, temperature))
        ofile.close()

    def _extractLast_PotEnergy(self,repl,cycle):
        replica = self.openmm_replicas[repl]
        (stateid, par) = replica.get_state()
        pot = replica.get_energy()
        pot_energy =  pot[0]
        temperature = par[0]
        if pot_energy == None:
            msg = "Error retrieving state for replica %d" % repl
            self._exit(msg)
        return (par, pot)

    def print_status(self):
        """
        Writes to BASENAME_stat.txt a text version of the status of the RE job

        It's fun to follow the progress in real time by doing:

        watch cat BASENAME_stat.txt
        """
        logfile = "%s_stat.txt" % self.basename
        ofile = self._openfile(logfile,"w")
        log = "Replica  State   Temperature Status  Cycle \n"
        for k in range(self.nreplicas):
            stateid = self.status[k]['stateid_current']
            log += "%6d   %5d  %s %5s  %5d \n" % (k, stateid, self.stateparams[stateid]['temperature'], self.status[k]['running_status'], self.status[k]['cycle_current'])
        log += "Running = %d\n" % self.running
        log += "Waiting = %d\n" % self.waiting

        ofile.write(log)
        ofile.close()

    def _getPot(self,repl,cycle):
        (par, pot) = self._extractLast_PotEnergy(repl,cycle)
        return pot

    def _getPar(self,repl):
        sid = self.status[repl]['stateid_current']
        tempt = float(self.stateparams[sid]['temperature'])
        kb = 0.0019872041
        beta = 1./(kb*tempt)
        return [beta]

    def _reduced_energy(self,par,pot):
        # par: list of parameters
        # pot: list of potentials
        # This is for temperature/binding potential beta*(U0+lambda*u)
        beta = par[0]
        epot = pot[0]
        return beta*epot

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

    rx = tempt_async_re_job(commandFile, options=None)

    rx.setupJob()

    rx.scheduleJobs()
