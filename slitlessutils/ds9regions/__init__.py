''' A package to write ds9-style regions file from within Python. '''


# import the attributes
from .info import __version__, __author__, __email__

# import the main routines
from .ds9regions import DS9Regions, ds9region


# import the regions:
from .annularregion import CircularAnnulus, EllipticalAnnulus, BoxAnnulus
from .astroregion import Compass, Ruler
from .compositeregion import Composite
from .enclosedregion import Circle, Ellipse, Box, Polygon
from .linearregion import Line, Vector, Projection, Segment
from .onedregion import Point, Text
# from .pandaregion import CircularPanda,EllipticalPanda,BoxPanda
