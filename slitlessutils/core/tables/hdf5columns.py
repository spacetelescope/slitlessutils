import numpy as np


class HDF5Columns(dict):
    ''' A primative class that will old the data for a table 


    has the numpy dtypes for each permissible column 

    ''' 

    DTYPES={'x': np.uint16,
            'y': np.uint16,
            'lam': np.uint16,
            'wav': np.float32,
            'xyg': np.uint64,
            'xyl': np.uint64,
            'val': np.float64}
 
