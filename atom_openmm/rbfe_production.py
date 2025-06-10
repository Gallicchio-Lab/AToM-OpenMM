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
import os

import openmm as mm
from openmm.app import *
from openmm import *
from openmm.unit import *
from datetime import datetime

from atom_openmm.openmm_async_re import openmm_job_RBFE
from atom_openmm.utils import set_directory
from pathlib import Path

if __name__ == '__main__':

    # Parse arguments:
    usage = "%prog <ConfigFile>"

    if len(sys.argv) != 2:
        print("Please specify ONE input file")
        sys.exit(1)

    commandFile = sys.argv[1]

    print("")
    print("========================================")
    print("AToM RBFE Asynchronous Replica Exchange ")
    print("========================================")
    print("")
    print("Started at: " + str(time.asctime()))
    print("Input file:", commandFile)
    print("")
    sys.stdout.flush()

    with set_directory(Path(commandFile).parent):
        rx = openmm_job_RBFE(os.path.basename(os.path.abspath(commandFile)), options=None)

        rx.setupJob()

        rx.scheduleJobs()
