import sys
import time

from sync.atm import openmm_job_AmberRBFE

if __name__ == '__main__':

    # Parse arguments:
    usage = "%prog <ConfigFile>"

    if len(sys.argv) != 2:
        print("Please specify ONE input file")
        sys.exit(1)

    config_file = sys.argv[1]

    print("")
    print("=======================================")
    print("AToM RBFE Synchronous Replica Exchange ")
    print("=======================================")
    print("")
    print("Started at: " + str(time.asctime()))
    print("Input file:", config_file)
    print("")
    sys.stdout.flush()

    rx = openmm_job_AmberRBFE(config_file)
    rx.setupJob()
    rx.scheduleJobs()
