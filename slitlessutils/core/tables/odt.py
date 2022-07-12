import numpy as np

from .hdf5table import HDF5Table

from ..utilities import indices

class ODT(HDF5Table):
    COLUMNS=('x','y','lam','val')
