import numpy as np


class HDF5Columns(dict):
    """
    Base class to hold the data for a given table

    Only contains a dictionary for the DTYPE for the possible columns

    """

    DTYPES = {'x': np.uint16,
              'y': np.uint16,
              'lam': np.uint16,
              'wav': np.float32,
              'xyg': np.uint64,
              'xyl': np.uint64,
              'val': np.float64}
