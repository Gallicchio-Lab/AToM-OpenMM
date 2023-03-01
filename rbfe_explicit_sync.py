import sys
from sync.atm import openmm_job_AmberRBFE

if __name__ == '__main__':
    assert len(sys.argv) == 2, "Specify ONE input file"

    rx = openmm_job_AmberRBFE(sys.argv[1])
    rx.setupJob()
    rx.scheduleJobs()
