import numpy as np
from scipy.constants import c

from ...logger import LOGGER


class Band:
    """
    Class to contain all sorts of filter bandpasses
    """

    def __init__(self, wave, tran, where=None, unit='angstrom', name=''):
        """
        Initializer

        Parameters
        ----------
        wave : `np.ndarray`
            the wavelength variable

        tran : `np.ndarray`
            the transmission array

        where : `np.ndarray` or None, optional
            A list of indices for the `wave` and `tran` to keep. if `None`
            then no indices will be subselected.  Default is `None`

        unit : str, optional
            the name of the units of the wavelength array. Default is 'angstrom'

        name : str, optional
            the name of the filter. Default is ''

        Notes
        -----
        Uses trapezoidal rule for numerical integration
        """

        if where is None:
            where = np.arange(len(wave), dtype=int)

        self.unit = unit.lower()
        self.wave = wave
        self.tran = tran
        self.name = name

        if self.unit in ('micron', 'microns', 'um'):
            self.wave *= 1e4
        elif self.unit in ('angstrom', 'angstroms', 'a', ''):
            pass
        else:
            LOGGER.warning(f'Unknown wavelength units {self.unit}. Assuming A')

        self.freq = (c/self.wave)*1e10

        # set ranges
        self.wmin = np.amin(self.wave[where])
        self.wmax = np.amax(self.wave[where])

        # compute the pivot wavelength
        # num=np.trapz(self.tran[where],x=self.wave[where])
        # den=np.trapz(self.tran[where]/self.wave[where]**2,x=self.wave[where])
        # self.photplam=np.sqrt(num/den)

        num = np.trapz(self.tran[where]*self.wave[where], x=self.wave[where])
        den = np.trapz(self.tran[where]/self.wave[where], x=self.wave[where])
        self.photplam = np.sqrt(num/den)

        # compute the max value
        self.tmax = np.amax(self.tran)

        # compute normalization
        self.fnunorm = np.trapz(self.tran/self.freq, x=self.freq)

    def __call__(self, l, left=0., right=0.):
        """
        Evaluate the bandpass by interpolation

        Parameters
        l : float, int, `np.ndarray`
            The wavelength to compute the transmission for

        left : float, optional
            Parameter passed to `np.interp()`, Default is 0.0

        right : float, optional
            Parameter passed to `np.interp()`, Default is 0.0

        Returns
        -------
        s : float or `np.ndarray'
            The interpolated sensitivity
        """

        s = np.interp(l, self.wave, self.tran, left=left, right=right)
        return s

    def __mul__(self, a):
        """
        Overide something to multiple a sensitivity curve by a scalar

        Parameters
        ----------
        a : float or int
           scalar to multiply by

        Returns
        -------
        a scaled version of the self
        """

        self.tran *= a
        self.tmax *= a
        return self

    def __rmul__(self, a):
        """
        see __mul__
        """
        return self.__mul__(a)

    def __imul__(self, a):
        """
        see __mul__
        """
        return self.__mul__(a)
