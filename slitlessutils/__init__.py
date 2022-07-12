# main imports
from .info import __code__
from .core import *
from . import examples

from .logger import initialize_logger
LOGGER=initialize_logger(__code__)
