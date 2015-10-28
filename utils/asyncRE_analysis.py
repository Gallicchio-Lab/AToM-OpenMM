# File Based Analysis Class for asynchronous Replica Exchang jobs
"""A module to analyze file-based asynchronous RE jobs
See documentation in doc/ directory.

Contributors:

Junchao Xia <junchao.xia@temple.edu>

"""

import os, sys, time
from configobj import ConfigObj
import random, math

def _exit(message):
    """Print and flush a message to stdout and then exit."""
    print message
    sys.stdout.flush()
    print 'exiting...'
    sys.exit(1)

class asyncRE_analysis:
    """
    Class to analyze results from asynchronous file-based RE calculations
    """

    def __init__(self, command_file, options):
        self.command_file = command_file
        self.jobname = os.path.splitext(os.path.basename(command_file))[0]
        self.keywords = ConfigObj(self.command_file)
        self._checkInput()
        self._printStatus()
        self._setupUWHAMTemplates()

    def _exit(self, message):
        _exit(message)

    def _printStatus(self):
        """Print a report of the input parameters."""
        print 'command_file =',self.command_file
        print 'jobname =',self.jobname
        for k,v in self.keywords.iteritems():
            print k,v


    def _checkInput(self):

        #Option for calculating binding free energy
        self.BindFreeEng=True
        if self.keywords.get('BIND_FREE_ENG') is None:
            self.BindFreeEng=False
        elif self.keywords.get('BIND_FREE_ENG').lower() == 'yes':
            self.BindFreeEng=True
        elif self.keywords.get('BIND_FREE_ENG').lower() == 'no':
            self.BindFreeEng=False
        else :
            self._exit("BIND_FREE_ENG option is not set right (yes or no).")
        if ( self.BindFreeEng) :
            if self.keywords.get('NBEGIN') is None:
                self._exit("The starting point (NBEGIN) for data needs to be specified")
            self.nbgn = int(self.keywords.get('NBEGIN'))
            if self.keywords.get('NEND') is None:
                self._exit("The end point (NEND) for data needs to be specified")
            self.nend = int(self.keywords.get('NEND'))
            if self.keywords.get('NFREQ') is None:
                self._exit("The frequency (NFREQ) for data needs to be specified")
            self.nfreq = int(self.keywords.get('NFREQ'))
        if self.keywords.get('CUMULATED') is None:
            self.cumulated=False
        elif self.keywords.get('CUMULATED').lower() == 'yes':
            self.cumulated=True
        elif self.keywords.get('CUMULATED').lower() == 'no':
            self.cumulated=False
        else :
            self._exit("CUMULATED option is not set right (yes or no).")
        if ( not self.cumulated) :
            if self.keywords.get('NDATA') is None:
                self._exit("The number of data point (NDATA) needs to be specified")
            self.ndata = int(self.keywords.get('NDATA'))


        #list of lambdas
        if self.keywords.get('LAMBDAS') is None:
            self._exit("LAMBDAS needs to be specified")
        self.lambdas = self.keywords.get('LAMBDAS').split(',')
        self.nlam = len(self.lambdas)

        #list of temperatures
        if self.keywords.get('TEMPERATURES') is None:
            self._exit("TEMPERATURES needs to be specified")
        self.temperatures = self.keywords.get('TEMPERATURES').split(',')
        self.ntemp = len(self.temperatures)
        #build parameters for the lambda/temperatures combined states
        self.nreplicas = self._buildBEDAMStates()

    def _buildBEDAMStates(self):
        self.stateparams = []
        for lambd in self.lambdas:
            for tempt in self.temperatures:
                st = {}
                st['lambda'] = lambd
                st['temperature'] = tempt
                self.stateparams.append(st)
        return len(self.stateparams)

    def _setupUWHAMTemplates(self):
        """ Setup templates for input files for R using UWHAM """

        self.uwham1D =  """

rm(list=ls())
npot.fcn <- function(x, lam) -bet*lam*x

# loads dataset
mydata = read.table("%s")
data(mydata)
lig.data <- mydata$V5
# sample size
N <- length(lig.data)
# lambda states
lam <- c(%s)
m <- length(lam)

# inverse temperature
bet <- 1.0/(0.001986209*%s)
# negative potential function


# state labels based on lambda values
# note that labels=1:m, not 0:(m-1)
state.labels <- factor(mydata$V4, labels=1:m)
# compute negative potential
neg.pot <- matrix(0, N,m)
for (j in 1:m)

neg.pot[,j] <- npot.fcn(x=lig.data, lam=lam[j])
# estimate free energies, note that size=NULL because label is given
require(UWHAM)
out <- uwham(label=state.labels, logQ=neg.pot, fisher=TRUE)
# free energies as a function of lambda, 0.36 kcal/mol is a standard
# state correction
ze <- matrix(out$ze, nrow=1, ncol=m)
-ze/bet
dg <- (-ze[,m]+ze[,1])/bet

# print out
# printf <- function(...)print(sprintf(...))
#
write(dg,file="%s", ncolumns = 1, append = FALSE, sep = " ")

# block bootstrap for free energies, note that proc.type="serial"
# for simulated tempering data.
# To save time for package checking, this is not run.
#out.boot <- uwham.boot(proc.type="serial", block.size=10, boot.size=100, label=state.labels, logQ=neg.pot)

#-out.boot$ze/bet
#sqrt(out.boot$ve)/bet

"""

        self.uwham2D =  """
library('trust')
library("UWHAM")

npot.fcn <- function(e0,ebind, bet, lam) -bet*(e0 + lam*ebind)

uwham.r <- function(label,logQ,ufactormax,ufactormin=1){
  n <- dim(logQ)[1]
  m <- dim(logQ)[2]
  iniz <- array(0,dim=m)
  uf <- ufactormax
  while(uf >= ufactormin & uf >= 1){
    mask <- seq(1,n,trunc(uf))
    out <- uwham(label=label.cross[mask], logQ=neg.pot[mask,],init=iniz)
    show(uf)
    iniz <- out$ze
    uf <- uf/2
  }
  out$mask <- mask
  out
}

data.t <- read.table("%s")
data.t$e0 <- data.t$V3 - data.t$V4 * data.t$V5
lam <- c(%s)
tempt <- c(%s)
bet <- 1.0/(0.001986209*tempt)
mtempt <- length(bet)
mlam <- length(lam)
m <- mlam*mtempt
N <- length(data.t$V1)

neg.pot <- matrix(0, N,m)
sid <- 1
# note the order of (be,te)
for (be in 1:mlam) {
     for (te in 1:mtempt) {
             neg.pot[,sid] <- npot.fcn(e0=data.t$e0,ebind=data.t$V5,bet[te],lam[be])
             sid <- sid + 1
    }
}
# note levels
label.tempt <- factor(data.t$V1, levels=tempt, labels=1:mtempt)
label.lam <- factor(data.t$V4, levels=lam, labels=1:mlam)
label.cross <- (as.numeric(label.lam)-1)*mtempt + as.numeric(label.tempt)
out <- uwham.r(label=label.cross, logQ=neg.pot,ufactormax=200,ufactormin=1)
ze <- matrix(out$ze, nrow=mtempt, ncol=mlam)
-ze/bet
dg <- (-ze[,mlam]+ze[,1])/bet
# print out
# printf <- function(...)print(sprintf(...))
# print
write(dg,file="%s", ncolumns = mtempt, append = FALSE, sep = " ")

"""
    def calculateBindFreeEng(self):
        """
calculate the binding free engies at different time from the time series of binding energies.
"""
        datafile = 'lbe_temp.dat'
        deltGfile = 'DeltG.dat'
        R_inpfile = "uwham_async.R"
        lambdas_str = ' '
        for il in range(0,self.nlam-1):
            lambdas_str +=  self.lambdas[il] + ','
        lambdas_str += self.lambdas[self.nlam-1]
        temps_str = ' '
        for it in range(0,self.ntemp-1):
            temps_str +=  self.temperatures[it] + ','
        temps_str += self.temperatures[self.ntemp-1]

        if (self.ntemp == 1) :
            uwham_input = self.uwham1D % (datafile,lambdas_str,temps_str,deltGfile)
        else:
            uwham_input = self.uwham2D % (datafile,lambdas_str,temps_str,deltGfile)
        f = open(R_inpfile, 'w')
        f.write(uwham_input)
        f.close()

        if (not self.cumulated):
            bfe_outf = open('bfe_conv_noc.dat', 'w')
        else:
            bfe_outf = open('bfe_conv.dat', 'w')

        for i in range(self.nbgn,self.nend+1):
            nhead = i*self.nfreq
            if (not self.cumulated):
                ndata = self.ndata
                ntail = self.ndata
            else:
                ndata = i*self.nfreq
                ntail = i*self.nfreq

            for ir in range(0,self.nreplicas):
                inpf = "r%d/lbe.dat" %ir
                outf = "r%d/lbe_temp.dat" %ir
                lbe_cmd = 'head -n ' + str(nhead) + ' ' + inpf + '| tail -n ' + str(ntail) + ' > ' + outf
                os.system(lbe_cmd)
                if (ir == 0 ):
                    tmp_cmd = 'cat ' + outf + '> lbe_temp.dat';
                else:
                    tmp_cmd = 'cat ' + outf + '>> lbe_temp.dat';
                os.system(tmp_cmd)
            uwham_cmd = 'R CMD BATCH ' + R_inpfile + '>& uwham_async.Rout'
            os.system(uwham_cmd)
            f = open(deltGfile ,"r")
            line = f.readline()
            f.close()
            #bfedata=line.split()
            #bfe_outf.write("%d" % int(ndata))
            #for j in range(0, len(bfedata)) :
            #    bfe_outf.write("%f" % float(bfedata[j]))
            #bfe_outf.write("\n")
            line = str(nhead) + '   ' + line
            bfe_outf.write(line)
        bfe_outf.close()



if __name__ == '__main__':

    # Parse arguments:
    usage = "%prog <ConfigFile>"

    if len(sys.argv) != 2:
        print "Please specify ONE input file"
        sys.exit(1)

    commandFile = sys.argv[1]

    print ""
    print "==============================================="
    print " Analyze Asynchronous Replica Exchange Results "
    print "==============================================="
    print ""
    print "Started at: " + str(time.asctime())
    print "Input file:", commandFile
    print ""
    sys.stdout.flush()

    async_analy = asyncRE_analysis(commandFile, options=None)

    if (async_analy.BindFreeEng) :
        async_analy.calculateBindFreeEng()
