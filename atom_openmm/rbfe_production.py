#! python

from __future__ import print_function
from __future__ import division
import sys
import time
import os
from multiprocessing import freeze_support

def rbfe_production(config_file=None, options=None):
    from atom_openmm.openmm_async_re import openmm_job_RBFE
    from atom_openmm.utils.AtomUtils import set_directory
    from pathlib import Path
    import os

    if config_file is None and options is None:
        config_file = sys.argv[1]

    assert config_file or options, "Invalid input. Specify a configuration file or options."

    print("")
    print("========================================")
    print("AToM RBFE Asynchronous Replica Exchange ")
    print("========================================")
    print("")
    print("Started at: " + str(time.asctime()))
    if config_file:
        print("Input file:", config_file)
    elif options:
        print("Options:", options)
    print("")
    sys.stdout.flush()

    work_dir = None
    
    if config_file:
        with set_directory(Path(config_file).parent):
            rx = openmm_job_RBFE(
                os.path.basename(os.path.abspath(config_file)), options=None
            )
            rx.setupJob()
            rx.scheduleJobs()
            work_dir = os.getcwd()
    elif options:
        if "WORKDIR" in options.keys():
            with set_directory(options["WORKDIR"]):
                rx = openmm_job_RBFE(config_file, options)
                rx.setupJob()
                rx.scheduleJobs()
                work_dir = os.getcwd()
        else:
            rx = openmm_job_RBFE(config_file, options)
            rx.setupJob()
            rx.scheduleJobs()
            work_dir = os.getcwd()

    return work_dir


if __name__ == "__main__":
    freeze_support()
    assert len(sys.argv) == 2, "Specify ONE input file"
    rbfe_production(config_file = sys.argv[1])
