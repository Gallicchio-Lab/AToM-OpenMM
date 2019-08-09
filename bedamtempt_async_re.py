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
    """
    def _doExchange_pair(self,repl_a,repl_b):
        #Performs exchange of lambdas for BEDAM replica exchange.

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
    """
    
    def _extractLast_lambda_BindingEnergy_TotalEnergy(self,repl,cycle):
        """
        Extracts binding energy from Impact output
        """
        output_file = "r%s/%s_%d.out" % (repl,self.basename,cycle)
        datai = self._getOpenMMData(output_file)
        nf = len(datai[0])
        nr = len(datai)
        # example of format (simple alchemical path):
        # <lambda, potential energy, binding energy>
        # 0.300000,-44753.762502,-130.875000
        #   0          1                 2
        #
        # example of format (quad bias alchemical path):
        # <lambda, gamma, b, w0, potential energy, binding energy>
        # 0.300000,0.1567,-1.345,12.345,-44753.762502,-130.875000
        #   0         1      2       3         4          5
        #format (ilogistic potential):
        # <lambda, lambda1, lambda2, alpha, u0, w0,  potential energy, binding energy>
        #   0         1       2       3     4    5       6                7
        #
        # [nr-1]: last record
        #
        if nf == 3:
            lmbd = datai[nr-1][0]
            binding_energy = datai[nr-1][2]
            pot_energy = datai[nr-1][1] #needs to be changed to total energy
            return (lmbd, binding_energy, pot_energy)
        elif nf == 6:
            lmbd =  datai[nr-1][0]
            gamma = datai[nr-1][1]
            b =     datai[nr-1][2]
            w0 =    datai[nr-1][3]
            parameters = [lmbd, gamma, b, w0]
            pot_energy = datai[nr-1][4] #needs to be changed to total energy
            binding_energy = datai[nr-1][5]
            return (parameters, binding_energy, pot_energy)
        elif nf == 8:
            lmbd =    datai[nr-1][0]
            lambda1 = datai[nr-1][1]
            lambda2 = datai[nr-1][2]
            alpha =   datai[nr-1][3]
            u0    =   datai[nr-1][4]
            w0    =   datai[nr-1][5]
            parameters = [lmbd, lambda1, lambda2, alpha, u0, w0]
            pot_energy = datai[nr-1][6] #needs to be changed to total energy
            binding_energy = datai[nr-1][7]
            return (parameters, binding_energy, pot_energy)

    def print_status(self):
        """
        Writes to BASENAME_stat.txt a text version of the status of the RE job

        It's fun to follow the progress in real time by doing:

        watch cat BASENAME_stat.txt
        """
        logfile = "%s_stat.txt" % self.basename
        ofile = self._openfile(logfile,"w")
        if self.gammas is not None:
            log = "Replica  State  Lambda Gamma Bcoeff W0coeff Temperature Status  Cycle \n"
            for k in range(self.nreplicas):
                stateid = self.status[k]['stateid_current']
                log += "%6d   %5d  %s %s %s %s %s %5s  %5d \n" % (k, stateid, self.stateparams[stateid]['lambda'], self.stateparams[stateid]['gamma'], self.stateparams[stateid]['bcoeff'], self.stateparams[stateid]['w0coeff'], self.stateparams[stateid]['temperature'], self.status[k]['running_status'], self.status[k]['cycle_current'])
        elif self.lambda1s is not None:
            log = "Replica  State  Lambda Lambda1 Lambda2 Alpha U0 W0coeff Temperature Status  Cycle \n"
            for k in range(self.nreplicas):
                stateid = self.status[k]['stateid_current']
                log += "%6d   %5d  %s %s %s %s %s %s %s %5s  %5d \n" % (k, stateid, self.stateparams[stateid]['lambda'], self.stateparams[stateid]['lambda1'], self.stateparams[stateid]['lambda2'], self.stateparams[stateid]['alpha'], self.stateparams[stateid]['u0'], self.stateparams[stateid]['w0coeff'], self.stateparams[stateid]['temperature'], self.status[k]['running_status'], self.status[k]['cycle_current'])
        else:
            log = "Replica  State  Lambda Temperature Status  Cycle \n"
            for k in range(self.nreplicas):
                stateid = self.status[k]['stateid_current']
                log += "%6d   %5d  %s %s %5s  %5d \n" % (k, stateid, self.stateparams[stateid]['lambda'], self.stateparams[stateid]['temperature'], self.status[k]['running_status'], self.status[k]['cycle_current'])
        log += "Running = %d\n" % self.running
        log += "Waiting = %d\n" % self.waiting

        ofile.write(log)
        ofile.close()

    def _getPot(self,repl,cycle):
        (parameters, u, etot) = self._extractLast_lambda_BindingEnergy_TotalEnergy(repl,cycle)
        if not isinstance(parameters, list): #lambda is the only alchemical parameter
            # removes lambda*u from etot to get e0. Note that this is the lambda from the
            # output file not the current lambda.
            lmb = parameters
            e0 = float(etot) - float(lmb)*float(u)
        elif len(parameters) == 4:
            lmb = float(parameters[0])
            gamma = float(parameters[1])
            b = float(parameters[2])
            w0 = float(parameters[3])
            uf = float(u)
            e0 = float(etot) - (0.5*gamma*uf*uf + b*uf + w0)
        elif len(parameters) == 6:
            lmb = float(parameters[0])
            lambda1 = float(parameters[1])
            lambda2 = float(parameters[2])
            alpha = float(parameters[3])
            u0 = float(parameters[4])
            w0 = float(parameters[5])
            uf = float(u)
            ee = 1.0 + math.exp(-alpha*(uf-u0))
            ebias = 0.0
            if alpha > 0:
                ebias = ((lambda2 - lambda1)/alpha) * math.log(ee)
            ebias += lambda2 * uf + w0
            e0 = float(etot) - ebias
        return (e0,uf)

    def _getPar(self,repl):
        sid = self.status[repl]['stateid_current']
        lmb = float(self.stateparams[sid]['lambda'])
        tempt = float(self.stateparams[sid]['temperature'])
        kb = 0.0019872041
        beta = 1./(kb*tempt)
        parameters = [beta,lmb]
        if 'gamma' in self.stateparams[sid].keys():
            gamma = float(self.stateparams[sid]['gamma'])
            parameters.append(gamma)
            b = float(self.stateparams[sid]['bcoeff'])
            parameters.append(b)
            w0 = float(self.stateparams[sid]['w0coeff'])
            parameters.append(w0)
        elif 'lambda1' in  self.stateparams[sid].keys():
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
        if len(par) == 2:
            # This is for temperature/binding potential beta*(U0+lambda*u)
            beta = par[0]
            lmb = par[1]
            e0 = pot[0]
            u = pot[1]
            return beta*(e0 + lmb*u)
        elif len(par) == 5:
            # This is for the quadbias potential
            # beta*(U0+0.5*gamma*u^2+b*u+w0)
            beta = par[0]
            lmb = par[1]
            gamma = par[2]
            b = par[3]
            w0 = par[4]
            e0 = pot[0]
            u = pot[1]
            return beta*(e0 + 0.5*gamma*u*u + b*u + w0)
        elif len(par) == 7:
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
        else:
            return None

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
