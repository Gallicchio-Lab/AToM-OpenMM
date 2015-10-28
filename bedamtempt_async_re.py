import sys
import time
import math
import random
import logging
from async_re import async_re
from bedam_async_re import bedam_async_re_job



class bedamtempt_async_re_job(bedam_async_re_job):
    def _setLogger(self):
        self.logger = logging.getLogger("async_re.bedamtempt_async_re")

    def _checkInput(self):
        async_re._checkInput(self)
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
        lambdas = self.keywords.get('LAMBDAS').split(',')
        #list of temperatures
        if self.keywords.get('TEMPERATURES') is None:
            self._exit("TEMPERATURES needs to be specified")
        temperatures = self.keywords.get('TEMPERATURES').split(',')
        #build parameters for the lambda/temperatures combined states
        self.nreplicas = self._buildBEDAMStates(lambdas,temperatures)
        #executive file's directory
        if self.keywords.get('JOB_TRANSPORT') is 'SSH':
            if self.keywords.get('EXEC_DIRECTORY') is None:
                self._exit("EXEC DIRECTORY needs to be specified")
        #added on 10.19.15
        self.implicitsolvent =  self.keywords.get('IMPLICITSOLVENT')
        self.totalsteps = self.keywords.get('TOTALSTEPS')
	self.jobname = self.keywords.get('ENGINE_INPUT_BASENAME')
	self.stepgap = self.keywords.get('STEPGAP')
	


    def _buildBEDAMStates(self,lambdas,temperatures):
        self.stateparams = []
        for lambd in lambdas:
            for tempt in temperatures:
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
        tbuffer = tbuffer.replace("@implicitsolvent@", implicitsolvent)
        tbuffer = tbuffer.replace("@totalsteps@", totalsteps)
	tbuffer = tbuffer.replace("@jobname@", jobname)
	tbuffer = tbuffer.replace("@stepgap@", stepgap)
        # write out
        ofile = self._openfile(inpfile, "w")
        ofile.write(tbuffer)
        ofile.close()

        # update the history status file
        ofile = self._openfile("r%d/state.history" % replica, "a")
        ofile.write("%d %d %s %s\n" % (cycle, stateid, lambd, temperature))
        ofile.close()

    def _doExchange_pair(self,repl_a,repl_b):
        """
        Performs exchange of lambdas for BEDAM replica exchange.
        """
        kb = 0.0019872041

        cycle_a = self.status[repl_a]['cycle_current']
        sid_a = self.status[repl_a]['stateid_current']
        lambda_a = float(self.stateparams[sid_a]['lambda'])
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
        (u_b,h_b) = self._extractLast_BindingEnergy_TotalEnergy(repl_b,cycle_b)
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
	#added on 10.21.15 check the exchange
	self.logger.info("dl = %f du = %f delta = %f", dl, du, delta)
	#added end on 10.21.15


        if self.keywords.get('VERBOSE') == "yes":
            self.logger.info("Pair Info")
            self.logger.info("%d %f %f %f %f", repl_a, lambda_a, u_a, beta_a, h_a)
            self.logger.info("%d %f %f %f %f", repl_b, lambda_b, u_b, beta_b, h_b)
            self.logger.info("dl = %f du = %f dh = %f delta = %f", dl, du, dh, delta)

        csi = random.random()
        if math.exp(-delta) > csi:
            status_func = lambda val: self.status[val]['stateid_current']
	    #added on 10.21.15 check the exchange
	    #self.logger.info("%s %s", status_func(repl_a), status_func(repl_b))
	    #added end on 10.21.15

            if self.keywords.get('VERBOSE') == "yes":
                self.logger.info("Accepted %f %f", math.exp(-delta), csi)
                self.logger.info("%s %s", status_func(repl_a), status_func(repl_b))
            self.status[repl_a]['stateid_current'] = sid_b
            self.status[repl_b]['stateid_current'] = sid_a
            if self.keywords.get('VERBOSE') == "yes":
                self.logger.info("%s %s", status_func(repl_a), status_func(repl_b))
        else:
            if self.keywords.get('VERBOSE') == "yes":
                self.logger.info("Rejected %f %f", math.exp(-delta), csi)

    def _extractLast_lambda_BindingEnergy_TotalEnergy(self,repl,cycle):
        """
        Extracts binding energy from Impact output
        """
        output_file = "r%s/%s_%d.out" % (repl,self.basename,cycle)
        datai = self._getOpenMMData(output_file)
        nf = len(datai[0])
        nr = len(datai)
        # [nr-1]: last record
        # [nf-2]: lambda (next to last item)
        # [nf-1]: binding energy (last item)
        #    [2]: total energy item (0 is step number and 1 is temperature)
        #
        # (lambda, binding energy, total energy)
	#check the binding energy in the last position 10.21.15
	#self.logger.info("binding energy in %s = %f ", output_file, datai[nr-1][nf-2])
	#end on 10.21.15
        return (datai[nr-1][nf-3],datai[nr-1][nf-2],datai[nr-1][nf-1]) #edit on 10.21.15

    def print_status(self):
        """
        Writes to BASENAME_stat.txt a text version of the status of the RE job

        It's fun to follow the progress in real time by doing:

        watch cat BASENAME_stat.txt
        """
        logfile = "%s_stat.txt" % self.basename
        ofile = self._openfile(logfile,"w")
        log = "Replica  State  Lambda Temperature Status  Cycle \n"
        for k in range(self.nreplicas):
            stateid = self.status[k]['stateid_current']
            log += "%6d   %5d  %s %s %5s  %5d \n" % (k, stateid, self.stateparams[stateid]['lambda'], self.stateparams[stateid]['temperature'], self.status[k]['running_status'], self.status[k]['cycle_current'])
        log += "Running = %d\n" % self.running
        log += "Waiting = %d\n" % self.waiting

        ofile.write(log)
        ofile.close()

    def _getPot(self,repl,cycle):
        (lmb, u, etot) = self._extractLast_lambda_BindingEnergy_TotalEnergy(repl,cycle)
        # removes lambda*u from etot to get e0. Note that this is the lambda from the
        # output file not the current lambda.
        e0 = float(etot) - float(lmb)*float(u)
        return (e0,float(u))

    def _getPar(self,repl):
        sid = self.status[repl]['stateid_current']
        lmb = float(self.stateparams[sid]['lambda'])
        tempt = float(self.stateparams[sid]['temperature'])
        kb = 0.0019872041
        beta = 1./(kb*tempt)
        return (beta,lmb)

    def _reduced_energy(self,par,pot):
        # par: list of parameters
        # pot: list of potentials
        # This is for temperature/binding potential beta*(U0+lambda*u)
        beta = par[0]
        lmb = par[1]
        e0 = pot[0]
        u = pot[1]
        return beta*(e0 + lmb*u)

if __name__ == '__main__':

    # Parse arguments:
    usage = "%prog <ConfigFile>"

    if len(sys.argv) != 2:
        print "Please specify ONE input file"
        sys.exit(1)

    commandFile = sys.argv[1]

    print ""
    print "===================================="
    print "BEDAM Asynchronous Replica Exchange "
    print "===================================="
    print ""
    print "Started at: " + str(time.asctime())
    print "Input file:", commandFile
    print ""
    sys.stdout.flush()

    rx = bedamtempt_async_re_job(commandFile, options=None)

    rx.setupJob()

    rx.scheduleJobs()
