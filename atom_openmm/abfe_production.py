#! python

from __future__ import print_function
from __future__ import division
import sys
import time


def abfe_production(config_file=None):
    from atom_openmm.openmm_async_re import openmm_job_ABFE
    from atom_openmm.utils.AtomUtils import set_directory
    from pathlib import Path
    import os

    if config_file is None:
        config_file = sys.argv[1]

    print("")
    print("=======================================")
    print("AToM ABFE Asynchronous Replica Exchange")
    print("=======================================")
    print("")
    print("Started at: " + str(time.asctime()))
    print("Input file:", config_file)
    print("")
    sys.stdout.flush()

    with set_directory(Path(config_file).parent):
        rx = openmm_job_ABFE(
            os.path.basename(os.path.abspath(config_file)), options=None
        )
        rx.setupJob()
        rx.scheduleJobs()


if __name__ == "__main__":
    assert len(sys.argv) == 2, "Specify ONE input file"

    abfe_production(sys.argv[1])
