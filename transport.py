"""A base class for all other transport modules (BOINC and SSH at the time
of writing) to inherit from.

The initial reason for this implementation is that async_re updates replica
statuses very inefficiently if the update process requires an external
connection to a database or something similar.  The initial version of the
boinc_transport module would open 3*N connections to the BOINC database every
runtime loop, where N is the number of replicas

Instead, we've added the Transport.poll() method which can be run to update the
status of all replicas at once.  For transport modules which don't need this
functionality, they will inherit an empty method
"""

import os
import logging, logging.config

class Transport(object):
    logging.config.fileConfig(os.path.join(os.path.dirname(__file__), "utils/logging.conf"))

    def __init__(self):
        pass

    def poll(self):
        return

