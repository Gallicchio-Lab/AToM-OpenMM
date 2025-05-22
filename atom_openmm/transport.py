import os
import logging, logging.config

class Transport(object):
    logging.config.fileConfig(os.path.join(os.path.dirname(__file__), "utils/logging.conf"))

    def __init__(self):
        pass

    def poll(self):
        return

