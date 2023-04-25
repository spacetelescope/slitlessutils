from dataclasses import dataclass

import numpy as np
from copy import copy

from ....logger import LOGGER
from ...utilities import headers


@dataclass
class Disperser:
    """
    Base dataclass to contain information on a dispersive element

    Parameters
    ----------
    name: str
        The name of the dispersive element

    blocking: str
        The name of the blocking element

    wave0: float or int
        The starting wavelength.

    wave1: float or int
        The ending wavelength

    units: str
        The units of the wavelength
    """

    name: str
    blocking: str
    wave0: float
    wave1: float
    units: str

    def __post_init__(self):
        """
        Method to check a few things

        Should never be directly called
        """

        wave0 = min(self.wave0, self.wave1)
        wave1 = max(self.wave0, self.wave1)

        self.wave0 = float(wave0)
        self.wave1 = float(wave1)

        if np.allclose(self.wave0, self.wave1):
            msg = 'wave0 and wave1 cannot be equal!: {self.wave0}, {self.wave1}'
            LOGGER.error(msg)

    def update_pars(self, specreg):
        """
        Method to update the parameters with locally set values for a
        spectral region, but any object with attributes will work

        Parameters
        ----------
        specreg : `su.core.sources.SpectralRegion`
           The spectral region to grab the parameters from.  Alternatively,
           this can be any object that has the correct attributes.

        Returns
        -------
        pars : `Disperser`
           A shallow copy of the class with updated values

        """

        # make the copy
        pars = copy(self)

        # if available, add the wavelength edges
        if hasattr(specreg, 'wave0'):
            pars.wave0 = specreg.wave0

        if hasattr(specreg, 'wave1'):
            pars.wave1 = specreg.wave1

        # based on the type, add the "delta" like parameters
        if isinstance(self, Linear) and hasattr(specreg, 'dwave'):
            pars.dwave = specreg.dwave
        elif isinstance(self, Geometric) and hasattr(specreg, 'scale'):
            pars.scale = specreg.scale

        # use this to validate the settings
        pars.__post_init__()

        return pars

    @property
    def trimmed_name(self):
        """
        Method to get a 'trimmed name'

        Returns
        -------
        name : str
            The trimmed name
        """
        print(self.__class__.__name__)
        if self.name[-1].isalpha():
            return self.name[:-1]
        else:
            return self.name

    def update_header(self, hdr):
        """
        Method to update a fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            The fits header to update
        """
        hdr.set('dispname', value=str(self.name), comment='Name of the disperser')
        hdr.set('blkname', value=str(self.blocking),
                comment='Name of the blocking filter')
        hdr.set('dispform', value=self.__class__.__name__,
                comment='Functional form of disperser')

        hdr.set('wave0', value=self.wave0, comment='Initial wavelength')
        hdr.set('wave1', value=self.wave1, comment='Final wavelength')
        hdr.set('wunit', value=self.units, comment='units on wavelength')

        headers.add_stanza(hdr, 'Wavelength Parameters', before='wave0')


@dataclass
class Linear(Disperser):
    """
    Dataclass to enable linear-dispersion spectral elements.

    inherits from `Disperser`

    Parameters
    ----------
    name: str
        The name of the dispersive element

    blocking: str
        The name of the blocking element

    wave0: float or int
        The starting wavelength.

    wave1: float or int
        The ending wavelength

    units: str
        The units of the wavelength

    dwave : float or int
        The step size in wavelength

    """

    dwave: float

    def __post_init__(self):
        """
        Method to check a few things

        Should never be directly called
        """

        assert (self.dwave > 0), 'Must have a positive wavelength sampling'
        super().__post_init__()
        self.dwave = float(self.dwave)

    def __len__(self):
        """
        Method to overload the len() function

        Returns
        -------
        nwav : int
            number of wavelength elements
        """

        nwav = int(np.ceil((self.wave1-self.wave0)/self.dwave)) + 1
        return nwav

    def __call__(self, lam, nsub=1):
        """
        Method to compute wavelengths from indices

        Parameters
        ----------
        lam : int or `np.ndarray` of ints
           The wavelength indices

        Returns
        -------
        wav : float or `np.ndarray` of floats
           The wavelengths
        """

        return self.wave0+lam*self.dwave/nsub

    def wavelengths(self, nsub=1):
        """
        Method to compute the center of the wavelength bins

        Parameters
        ----------
        nsub : int, optional
            The subsampling frequency.  Default is 1.

        Returns
        -------
        wav : `np.ndarray`
            The center of the wavelengths bins
        """

        wav = np.arange(self.wave0, self.wave1+self.dwave/nsub,
                        self.dwave/nsub, dtype=float)
        return wav

    def indices(self, wav):
        """
        Method to convert from floating-point wavelengths to indices

        Parameters
        ----------
        wav : float, int, `np.ndarray`
            The wavelengths to compute

        Returns
        -------
        ind : int, `np.ndarray`
            The indices (will be integer dtype)
        """

        ind = np.floor((wav-self.wave0)/self.dwave).astype(int)
        return ind

    def limits(self, nsub=1):
        """
        Method to compute the wavelength limits

        Parameters
        ----------
        nsub : int, optional
            The subsampling frequency.  Default is 1.

        Returns
        -------
        lim : `np.ndarray`
            A float array that gives the edges of the bins
        """
        dw2 = self.dwave/2.
        dwn = self.dwave/float(nsub)

        lim = np.arange(self.wave0-dw2, self.wave1+dw2+dwn, dwn)
        return lim

    def update_header(self, hdr):
        """
        Method to update a fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            The fits header to update
        """

        super().update_header(hdr)
        hdr.set('dwave', value=self.dwave, after='wave1',
                comment='Sampling wavelength')


@dataclass
class Geometric(Disperser):
    """
    Dataclass to enable geometric-dispersion spectral elements.

    inherits from `Disperser`

    Parameters
    ----------
    name: str
        The name of the dispersive element

    blocking: str
        The name of the blocking element

    wave0: float or int
        The starting wavelength.

    wave1: float or int
        The ending wavelength

    units: str
        The units of the wavelength

    scale : float or int
        The multiplicative scale

    """

    scale: float

    def __post_init__(self):
        """
        Method to check a few things

        Should never be directly called
        """

        assert (self.scale > 1), 'Must have an increasing scale factor'
        super().__post_init__()
        self.scale = float(self.scale)

    def __len__(self):
        """
        Method to overload the len() function

        Returns
        -------
        nwav : int
            number of wavelength elements
        """
        pass

    def __call__(self, lam, nsub=1):
        pass

    def wavelengths(self, nsub=1):
        """
        Method to compute the center of the wavelength bins

        Parameters
        ----------
        nsub : int, optional
            The subsampling frequency.  Default is 1.

        Returns
        -------
        wav : `np.ndarray`
            The center of the wavelengths bins
        """
        pass

    def indices(self, wav):
        """
        Method to convert from floating-point wavelengths to indices

        Parameters
        ----------
        wav : float, int, `np.ndarray`
            The wavelengths to compute

        Returns
        -------
        ind : int, `np.ndarray`
            The indices (will be integer dtype)
        """
        pass

    def limits(self, nsub=1):
        """
        Method to compute the wavelength limits

        Parameters
        ----------
        nsub : int, optional
            The subsampling frequency.  Default is 1.

        Returns
        -------
        lim : `np.ndarray`
            A float array that gives the edges of the bins
        """

        pass

    def update_header(self, hdr):
        """
        Method to update a fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            The fits header to update
        """

        super().update_header(hdr)
        hdr.set('scale', value=self.scale, after='wave1',
                comment='logarithmic sampling')


def load_disperser(name, blocking, **kwargs):
    disptype = kwargs['disptype'].lower()
    if disptype == 'grism':
        disp = Linear(name, blocking, kwargs['wave0'], kwargs['wave1'],
                      kwargs['units'], kwargs['dwave'])
    elif disptype == 'prism':
        disp = Geometric(name, blocking, kwargs['wave0'], kwargs['wave1'],
                         kwargs['units'], kwargs['scale'])
    else:
        raise NotImplementedError(f'Disperser type {disptype} is unknown.')

    return disp
