from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("atom-openmm")
except PackageNotFoundError:
    pass
