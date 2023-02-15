import sys
import time

from openmm_sync_re import openmm_job_AmberRBFE

if __name__ == '__main__':

    # Parse arguments:
    usage = "%prog <ConfigFile>"

    if len(sys.argv) != 2:
        print("Please specify ONE input file")
        sys.exit(1)

    commandFile = sys.argv[1]

    print("")
    print("=======================================")
    print("AToM RBFE Synchronous Replica Exchange ")
    print("=======================================")
    print("")
    print("Started at: " + str(time.asctime()))
    print("Input file:", commandFile)
    print("")
    sys.stdout.flush()

    rx = openmm_job_AmberRBFE(commandFile, options=None)
    rx.setupJob()
    rx.scheduleJobs()
