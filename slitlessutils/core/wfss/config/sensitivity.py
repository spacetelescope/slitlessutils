import os
import numpy as np

from astropy.io import fits

from ...photometry import Band


class Sensitivity(Band):
    """
    Class that holds the WFSS sensitivity curve.

    inherits from `su.core.photometry.Band`
    """

    def __init__(self, sensfile, senslimit=1e10):
        """
        Initialize the `Sensitivity` object.

        Parameters
        ----------
        sensfile : str
            The full path to the sensitivity file, which should be a fits
            table.

        senslimit : float
            The limit, below which, the data are ignored for the purposes
            of computing the average values.  Default is 1e10.

        """
        self.sensfile = sensfile
        self.senslimit = senslimit

        # read the file
        if os.path.exists(self.sensfile):
            data, header = fits.getdata(self.sensfile, 1, header=True)

            g = np.where(data['SENSITIVITY'] > self.senslimit)
            Band.__init__(self, data['WAVELENGTH'], data['SENSITIVITY'],
                          where=g, unit=header.get('TUNIT1', ''))
        else:
            Band.__init__(self, np.array([0, np.inf]), np.array([1., 1.]))

    def __str__(self):
        return f'Sensitivity curve: {self.sensfile}'
