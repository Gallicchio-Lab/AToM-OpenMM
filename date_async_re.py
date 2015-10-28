import os, sys, time
from async_re import async_re

"""
Adapted from:
https://github.com/saga-project/asyncre-bigjob
"""

class date_job(async_re):

    def _launchReplica(self,replica,cycle):
        """
        Issues a command to launch /bin/date
        """
        job_info = {
            "executable": "dodate",
            "input_file": "",
            "output_file": "sj-stdout-"+str(replica)+"-"+str(cycle)+".txt",
            "error_file": "sj-stderr-"+str(replica)+"-"+str(cycle)+".txt",
            "working_directory":os.getcwd()+"/r"+str(replica),
        }

        if self.keywords.get('VERBOSE') == "yes":
            print "Launching %s in directory %s cycle %d" % ("/bin/date",os.getcwd()+"/r"+str(replica),cycle)

        status = self.transport.launchJob(replica, job_info)

        return status


class date_async_re_job(date_job,async_re):

    def _checkInput(self):
        async_re._checkInput(self)
        #make sure DATE is wanted
        if self.keywords.get('RE_TYPE') != 'DATE':
            self._exit("RE_TYPE is not DATE")
        #number of replicas
        if self.keywords.get('NREPLICAS') is None:
            self._exit("NREPLICAS needs to be specified")
        self.nreplicas = int(self.keywords.get('NREPLICAS'))

    def _hasCompleted(self,replica,cycle):
        return True

    def _buildInpFile(self, replica):
        pass

    def _doExchange_pair(self,repl_a,repl_b):
        pass

    def _computeSwapMatrix(self, replicas, states):
        U = [[ 0. for j in range(self.nreplicas)]
             for i in range(self.nreplicas)]
        return U

if __name__ == '__main__':

    # Parse arguments:
    usage = "%prog <ConfigFile>"

    if len(sys.argv) != 2:
        print "Please specify ONE input file"
        sys.exit(1)

    commandFile = sys.argv[1]

    print ""
    print "===================================="
    print "DATE Asynchronous Replica Exchange "
    print "===================================="
    print ""
    print "Started at: " + str(time.asctime())
    print "Input file:", commandFile
    print ""
    sys.stdout.flush()

    rx = date_async_re_job(commandFile, options=None)

    rx.setupJob()

    rx.scheduleJobs()
