import os
import sys
import re
import time
import math
import random
import logging
from async_re import async_re
from impact_async_re import impact_job

class tempt_async_re_job(impact_job):

    def _setLogger(self):
        self.logger = logging.getLogger('async_re.tempt_async_re')


    def _checkInput(self):
        async_re._checkInput(self)
        #make sure TEMPERATURE is wanted
        if self.keywords.get('RE_TYPE') != 'TEMPT':
            self._exit("RE_TYPE is not TEMPT")
        #we run with IMPACT
        if self.keywords.get('ENGINE') != 'IMPACT':
            self._exit("ENGINE is not IMPACT")
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
        basename = self.basename
        stateid = self.status[replica]['stateid_current']
        cycle = self.status[replica]['cycle_current']

        template = "%s.inp" % basename
        inpfile = "r%d/%s_%d.inp" % (replica, basename, cycle)

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
        tbuffer = tbuffer.replace("@cycle@",str(cycle))
        # write out
        ofile = self._openfile(inpfile, "w")
        ofile.write(tbuffer)
        ofile.close()

        # update the history status file
        ofile = self._openfile("r%d/state.history" % replica, "a")
        ofile.write("%d %d %s\n" % (cycle, stateid, temperature))
        ofile.close()

    def _doExchange_pair(self,repl_a,repl_b):
        """
        Performs exchange of lambdas for BEDAM replica exchange.
        """
        kb = 0.0019872041

        cycle_a = self.status[repl_a]['cycle_current']
        sid_a = self.status[repl_a]['stateid_current']
        temperature_a = float(self.stateparams[sid_a]['temperature'])
        # u_a: binding energy of replica a
        # h_a: total energy of replica a (includes kinetic energy as we are not
        #      doing velocity rescaling here)
        (u_a,h_a) = self._extractLast_BindingEnergy_TotalEnergy(repl_a,cycle_a)
        (u_a,h_a) = (float(u_a),float(h_a))

        cycle_b = self.status[repl_b]['cycle_current']
        sid_b = self.status[repl_b]['stateid_current']
        lambda_b = float(self.stateparams[sid_b]['lambda'])
        temperature_b = float(self.stateparams[sid_b]['temperature'])
        (u_b,h_b) = self._extractLast_TotalEnergy(repl_b,cycle_b)
        (u_b,h_b) = (float(u_b),float(h_b))

        # Acceptance criterion is based on exp(-Delta) where
        #  Delta = -(beta_b - beta_a)*[H_b-H_a] - (lmbd_b - lmbd_a)[beta_a*u_b-beta_b*u_a]
        # To derive this start from the Boltzmann weight exp[-F(x|lambda,beta)]
        # where F(x|lambda,beta) = beta*[H_0(x)+lambda*u(x)], set up the usual
        # Metropolis exchange rules noticing that:
        #  H_b(x_a) = H_a(x_a) + (lmbd_b - lmbd_a)*u(x_a)

        beta_a = 1./(kb*temperature_a)
        beta_b = 1./(kb*temperature_b)
        dl = lambda_b - lambda_a
        du = beta_a*u_b - beta_b*u_a
        db = beta_b - beta_a
        dh = h_b - h_a
        delta = -dl*du - db*dh

        if self.keywords.get('VERBOSE') == "yes":
            self.logger.info("Pair Info")
            self.logger.info("%d %f %f %f %f", repl_a, lambda_a, u_a, beta_a, h_a)
            self.logger.info("%d %f %f %f %f", repl_b, lambda_b, u_b, beta_b, h_b)
            self.logger.info("dl = %f du = %f dh = %f delta = %f", dl, du, dh, delta)

        csi = random.random()
        if math.exp(-delta) > csi:
            status_func = lambda val: self.status[val]['stateid_current']
            if self.keywords.get('VERBOSE') == "yes":
                self.logger.info("Accepted %f %f", math.exp(-delta), csi)
                self.logger.info("%s %s", status_func(repl_a), status_func(repl_b))
            self.status[repl_a]['stateid_current'] = sid_b
            self.status[repl_b]['stateid_current'] = sid_a
            if self.keywords.get('VERBOSE') == "yes":
                self.logger.info("%s %s", status_func(repl_a), status_func(repl_b))
        else:
            if self.keywords.get('VERBOSE') == "yes":
                self.logger.info("Accepted %f %f", math.exp(-delta), csi)

    def _getImpactData(self, file):
        """
        Override because with T-RE there are only 2 lines of data printed in the
        output file.

        Reads all of the Impact simulation data values temperature, energies,
        etc.  at each time step and puts into a big table
        """
        if not os.path.exists(file):
            msg = 'File does not exist: %s' % file
            self._exit(msg)
        step_line = re.compile("^ Step number:")
        number_line = re.compile("(\s+-*\d\.\d+E[\+-]\d+\s*)+")
        nsamples = 0
        data = []
        f = self._openfile(file ,"r")
        line = f.readline()
        while line:
            # fast forward until we get to the line:
            # "Step number: ... "
            while line and not re.match(step_line, line):
                line = f.readline()
            # read the step number
            if re.match(step_line, line):
                words = line.split()
                step = words[2]
                #now read up to 2 lines of numbers
                datablock = [int(step)]
                ln = 0
                while ln < 2:
                    line = f.readline()
                    if not line:
                        msg = "Unexpected end of file"
                        self._exit(msg)
                    if re.match(number_line, line):
                        for word in line.split():
                            datablock.append(float(word))
                        ln += 1
                data.append(datablock)
            line = f.readline()
        f.close()
        return data

    def _extractLast_TotalEnergy(self,repl,cycle):
        """
        Extracts binding energy from Impact output
        """
        output_file = "r%s/%s_%d.out" % (repl,self.basename,cycle)
        datai = self._getImpactData(output_file)
        nf = len(datai[0])
        nr = len(datai)
        # [nr-1]: last record
        #    [2]: total energy item (0 is step number and 1 is temperature)
        #
        # (total energy)
        return datai[nr-1][2]

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
        etot = self._extractLast_TotalEnergy(repl,cycle)
        return [etot]

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
        e0 = pot[0]
        return beta*e0

if __name__ == '__main__':

    # Parse arguments:
    usage = "%prog <ConfigFile>"

    if len(sys.argv) != 2:
        print "Please specify ONE input file"
        sys.exit(1)

    commandFile = sys.argv[1]

    print ""
    print "===================================="
    print "Temperature Asynchronous Replica Exchange "
    print "===================================="
    print ""
    print "Started at: " + str(time.asctime())
    print "Input file:", commandFile
    print ""
    sys.stdout.flush()

    rx = tempt_async_re_job(commandFile, options=None)

    rx.setupJob()

    rx.scheduleJobs()
