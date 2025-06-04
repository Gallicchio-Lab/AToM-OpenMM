from importlib.metadata import version, PackageNotFoundError
import logging.config
import os

__curr_dir = os.path.dirname(os.path.abspath(__file__))

try:
    __version__ = version("atom-openmm")
except PackageNotFoundError:
    pass


try:
    logging.config.fileConfig(
        os.path.join(__curr_dir, "utils", "logging.conf"),
        disable_existing_loggers=False,
    )
except Exception as e:
    print(f"AToM-OpenMM: Logging setup failed: {e}")
