# main imports
from .info import __version__   # noqa: F401
from .core import *  # noqa: F401, F403
from .ds9regions import *  # noqa: F401, F403
from . import examples  # noqa: F401

from .logger import start_logging, end_logging, initialize_logger  # noqa: F401
LOGGER = initialize_logger()
