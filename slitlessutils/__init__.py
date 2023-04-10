# main imports
from .info import __version__#__code__
from .core import *
from .ds9regions import *
from . import examples

from .logger import start_logging,end_logging,initialize_logger
LOGGER=initialize_logger()
