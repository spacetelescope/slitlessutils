# main imports
from . import examples  # noqa: F401
from .core import *  # noqa: F401, F403
from .ds9regions import *  # noqa: F401, F403

try:
    from .version import version as __version__
except ImportError:
    __version__ = ''
