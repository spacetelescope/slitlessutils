
try:
    from .version import version as __version__
except ImportError:
    __version__ = ''

from importlib.metadata import metadata

__code__ = 'slitlessutils'
package_info = metadata(__code__)

__author__, __email__ = package_info['Author-email'].split(' <')
__email__ = __email__.replace('>', '')
