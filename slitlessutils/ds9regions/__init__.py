''' A package to write ds9-style regions file from within Python. '''


# import the attributes
from .info import __version__, __author__, __email__  # noqa: F401

# import the main routines
from .ds9regions import DS9Regions, ds9region  # noqa: F401

# import the regions:
from .annularregion import CircularAnnulus, EllipticalAnnulus, BoxAnnulus  # noqa: F401
from .astroregion import Compass, Ruler  # noqa: F401
from .compositeregion import Composite  # noqa: F401
from .enclosedregion import Circle, Ellipse, Box, Polygon  # noqa: F401
from .linearregion import Line, Vector, Projection, Segment  # noqa: F401
from .onedregion import Point, Text  # noqa: F401
# from .pandaregion import CircularPanda,EllipticalPanda,BoxPanda  # noqa: F401
