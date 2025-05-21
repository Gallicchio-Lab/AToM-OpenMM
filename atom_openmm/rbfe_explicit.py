#! python

from __future__ import print_function
from __future__ import division
import sys
import time
import math
import random
import logging
import signal
import shutil
import random

import openmm as mm
from openmm.app import *
from openmm import *
from openmm.unit import *
from datetime import datetime

from atom_openmm.openmm_async_re import openmm_job_RBFE

if __name__ == '__main__':

    # Parse arguments:
    usage = "%prog <ConfigFile>"

    if len(sys.argv) != 2:
        print("Please specify ONE input file")
        sys.exit(1)

    commandFile = sys.argv[1]

    print("")
    print("DEPRECATION WARNING:")
    print("The rbfe_explicit.py script is no longer supported")
    print("and will be removed in future releases.")
    print("Use rbfe_production.py instead.")
    print("")

    print("")
    print("========================================")
    print("AToM RBFE Asynchronous Replica Exchange ")
    print("========================================")
    print("")
    print("Started at: " + str(time.asctime()))
    print("Input file:", commandFile)
    print("")
    sys.stdout.flush()

    rx = openmm_job_RBFE(commandFile, options=None)

    rx.setupJob()

    rx.scheduleJobs()
