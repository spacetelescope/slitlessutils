''' A package to write ds9-style regions file from within Python. '''


# import the regions:
from .annularregion import (BoxAnnulus, CircularAnnulus,  # noqa: F401
                            EllipticalAnnulus)
from .astroregion import Compass, Ruler  # noqa: F401
from .compositeregion import Composite  # noqa: F401
# import the main routines
from .ds9regions import DS9Regions, ds9region  # noqa: F401
from .enclosedregion import Box, Circle, Ellipse, Polygon  # noqa: F401
# import the attributes
from .info import __author__, __email__, __version__  # noqa: F401
from .linearregion import Line, Projection, Segment, Vector  # noqa: F401
from .onedregion import Point, Text  # noqa: F401

# from .pandaregion import CircularPanda,EllipticalPanda,BoxPanda  # noqa: F401
