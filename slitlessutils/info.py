__code__ = 'slitlessutils'

from ._version import version as __version__  # noqa: F401
from importlib.metadata import metadata

d = metadata(__code__)

__author__ = d['Author']
__email__ = d['Email']
__version__ = d['Version']
