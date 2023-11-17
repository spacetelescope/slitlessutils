import os
from dataclasses import dataclass
from datetime import datetime

import numpy as np
import pysiaf
import yaml
from astropy.io import fits
from pypolyclip import clip_multi

from ....config import Config
from ....logger import LOGGER
from ...utilities import headers
from .disperser import Disperser, load_disperser
from .wfssconfig import WFSSConfig

# MJD is defined as number of days since midnight on November 17, 1858
MJDREF = datetime(1858, 11, 17, 0, 0, 0)

# conversion table between strings and datatypes
DTYPES = {'np.float16': np.float16,
          'np.float32': np.float32,
          'np.float64': np.float64,
          'np.int8': np.int8,
          'np.int16': np.int16,
          'np.int32': np.int32,
          'np.int64': np.int64,
          'np.uint8': np.uint8,
          'np.uint16': np.uint16,
          'np.uint32': np.uint32,
          'np.uint64': np.uint64,
          'int': int,
          'float': float}


@dataclass
class SIAF:
    """
    Dataclass to contain information for the pointing of a detector or
    instrument.

    One should probably never directly use this class.

    Parameters
    ----------
    v2 : float
        The V2 coordinate (in deg)

    v3 : float
        The V3 coordinate (in deg)

    v3y : float
        The rotation of the V3 axis (in deg)
    """
    v2: float
    v3: float
    v3y: float

    def update_header(self, hdr, reference=False):
        """
        Method to update a fits header

        Parameters
        ----------
        hdr : `astropy.fits.io.Header`
            Fits header to update

        reference : bool, optional
            Flag if this is for a reference SIAF (ie. one that describes
            the instrument pointing, as opposed to a detector's pointing)
        """

        ref = 'Reference' if reference else ''
        hdr.set('V2', value=self.v2, comment='V2 (arcsec)')
        hdr.set('V3', value=self.v3, comment='V3 (arcsec)')
        hdr.set('V3Y', value=self.v3y, comment='V3Y (deg)')
        headers.add_stanza(hdr, f'{ref} SIAF Info', before='V2')

    def compute_pav3(self, orientat):
        """
        Method to compute the PA_V3

        Parameters
        ----------
        orientat : float
            The orientation.  Angle East of North in deg

        Returns
        -------
        pav3 : float
            The position angle of the v3 axis
        """

        pav3 = (orientat - self.v3y) % 360

        return pav3

    def compute_crvals(self, ra, dec, orientat, refsiaf):
        """
        Method to compute the CRVALs given a pointing and a reference SIAF
        for the instrument

        Parameters
        ----------
        ra : float
            The RA of the field center (in deg)

        dec : float
            The Dec of the field center (in deg)

        orientat : float
            The position angle (East of North) of the field (in deg)

        refsiaf : `SIAF`
            An instance of the `SIAF` class that describes the instrument

        Returns
        -------
        crvals : 2-tuple
            A tuple of the two CRVALs (in deg)

        """

        # do a quick test, if the refsiaf == self, then
        # do not need to do anything
        if refsiaf == self:
            # if here, then the reference data is where we're at.  so,
            # let's just short-circute and use the ra/dec as given
            crvals = (ra, dec)
        else:
            # if here, then have to compute the CRVALS, so do this
            # based on Bryan Hilbert's email from Apr 15, 2020, I should
            # transform orientat to a local_roll, but that is really
            # complicated.  But the first order is:
            # local_roll=-orientat-self.v3y
            # is essentially that.  Higher order will require more coding.
            # local_roll=-orientat-self.v3y
            local_roll = orientat - self.v3y
            A = pysiaf.rotations.attitude(refsiaf.v2, refsiaf.v3,
                                          ra, dec, local_roll)

            # compute the new positions
            crvals = pysiaf.utils.rotations.pointing(A, self.v2, self.v3)
        return crvals

    def __eq__(self, siaf):
        if isinstance(siaf, type(self)):
            return np.isclose(self.v2, siaf.v2) and \
                np.isclose(self.v3, siaf.v3) and \
                np.isclose(self.v3y, siaf.v3y)
        else:
            return False


@dataclass
class Extension:
    """
    A dataclass to hold information for a fits-file extension

    Parameters
    ----------
    extname : str
        The name of the extension (will be EXTNAME in the fits header)

    extver : int
        The version of the extension (will be EXTVER in the fits header)

    dtype : type
        A valid python type indicator (used to recast data as necessary).
        If an invalid type is set, then `float` will be assumed.

    """

    extname: str
    extver: int
    dtype: type

    def __post_init__(self):
        """
        A validation method.

        Should never be directly called
        """
        self.dtype = DTYPES.get(self.dtype, 'float')

        # pep doesn't like this:  so deprecating this...
        # try:
        #     self.dtype = eval(self.dtype)
        # except BaseException:
        #     self.dtype = float

    @property
    def extension(self):
        """
        The tuple of the extension label

        Returns
        -------
        ext : 2-tuple
            The (extname,extver) to pass to `astropy.io.fits` codes
        """
        ext = (self.extname, self.extver)
        return ext

    def update_header(self, hdr):
        """
        Method to update a fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            The fits header to update
        """

        hdr['EXTNAME'] = (self.extname, "Name of extension")
        hdr['EXTVER'] = (self.extver, 'Extension version')

    def retype(self, dat):
        """
        Method to retype a datum according to the set dtype

        Parameters
        ----------
        dat : `np.ndarray`
            The data to retype

        Returns
        -------
        The retyped data
        """
        return dat.astype(self.dtype)

    def readfits(self, filename, header=False):
        """
        Method to load data from a fits file

        The extension according to this `Extension` object is read.

        Parameters
        ----------
        filename : str
            The filename to read

        header : bool, optional
            Flag to also read and return the extension header.
            Default is False

        Returns
        -------
        The image and/or the the header as requested by the `header` flag.

        """

        # with warnings.catch_warnings():
        #    warnings.simplefilter('ignore',FITSFixedWarning)
        a = fits.getdata(filename, extname=self.extname,
                         extver=self.extver, header=header)
        return a

    def headfits(self, filename):
        """
        Method to load header from a fits file

        The extension according to this `Extension` object is read.

        Parameters
        ----------
        filename : str
            The filename to read

        Returns
        -------
        h : `astropy.io.fits.Header`
            The header

        """
        # with warnings.catch_warnings():
        #    warnings.simplefilter('ignore',FITSFixedWarning)
        h = fits.getheader(filename, extname=self.extname, extver=self.extver)
        return h


@dataclass
class Noise:
    """
    A Dataclass to handle all things regarding noise calculations

    Parameters
    ----------
    dark : float
        The dark rate in e-/s

    read : float
        The read noise in e-

    time : float
        The notional exposure time in s

    back : float
        The approximate background rate in e-/s

    Notes
    -----
    This module is likely to change in the future as new instruments and
    more complex modeling becomes necesssary or available.
    """

    dark: float
    read: float
    time: float
    back: float

    def __call__(self, sci):
        """
        Method to add noise to an image

        Parameters
        ----------
        sci : float, int, `np.ndarray`
            The science data to add noise to (in e-/s)

        Returns
        -------
        new : float, int, `np.ndarray`
            The updated science image that has noise (in e-/s)

        unc : float, int, `np.ndarray`
            The uncertainty image (in e-/s)

        Notes
        -----
        The output images will have same shape/dtype as the input science
        image.
        """

        pvars = sci * self.time
        pvarb = (self.dark + self.back) * self.time
        new = np.random.poisson(lam=pvars + pvarb, size=sci.shape) - pvarb

        gvar = self.read**2
        new = np.random.normal(loc=new, scale=np.sqrt(gvar)) / self.time
        unc = np.sqrt(pvars + pvarb + gvar) / self.time

        return new, unc

    def update_header(self, hdr):
        """
        A method to update a fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            The fits header to update
        """
        hdr['NOISE'] = (True, 'was noise added?')
        hdr['RDNOISE'] = (self.read, 'assumed read noise in e-')
        hdr['DARKRATE'] = (self.dark, 'assumed dark rate in e-/2')
        headers.add_stanza(hdr, 'Noise Properties', before='NOISE')


class DetectorConfig:
    """
    Class to configure a WFSS detector
    """

    def __init__(self, name, data, disperser, path, exptime=1000., background=0.,
                 **kwargs):
        """
        Initializer method

        Parameters
        ----------
        name : str
            The name of the detector

        data : dict
            A dictionary that contains information from the yaml config file

        disperse : `Disperser`
            A `Disperser` object to configure

        path : str
            A path to the reference files

        exptime : float, optional
            The notional exposure time (in sec) for noise calculations.
            Default is 1000.

        background : float, optional
            The background level in (e-/s) for noise calculations. Default
            is 0.0

        kwargs : dict, optional
            optional keywords that get passed to `WFSSConfig`

        """

        # record a handful of things
        self.name = name
        self.naxis = np.array(data['naxis'], dtype=int)
        self.crpix = np.array(data['crpix'], dtype=float)
        self.scale = np.array(data['scale'], dtype=float)
        self.noise = Noise(*data['noise'].values(), exptime, background)
        self.header = data.get('header', {})
        self.bitmask = data.get('bitmask')

        self.extensions = {k: Extension(*v.values()) for k, v in
                           data['extensions'].items()}

        self.siaf = SIAF(*data['siaf'].values())

        if data['blockers'][disperser.blocking]['sip']:
            a = self.getSipMatrix(data['blockers'][disperser.blocking]['sip']['a'])
            b = self.getSipMatrix(data['blockers'][disperser.blocking]['sip']['b'])

            self.sip = (a, b, None, None, self.crpix)
        else:
            self.sip = ()

        conffile = data['blockers'][disperser.blocking]['config'][disperser.name]
        conffile = os.path.join(Config().refpath, 'instruments', path, conffile)
        self.config = WFSSConfig(conffile, **kwargs)

    @staticmethod
    def getSipMatrix(data):
        """
        Staticmethod to extract SIP data from a yaml-file dict

        Parameters
        ----------
        data : dict
            A dict from the yaml file

        Returns
        -------
        sip : `np.ndarray`
            The SIP coefficients as a 2d array.

        """

        i = []
        j = []
        v = []
        for key, val in data.items():
            ii, jj = key.split(',')
            i.append(int(ii))
            j.append(int(jj))
            v.append(float(val))
        n = np.maximum(np.amax(i), np.amax(j)) + 1
        sip = np.zeros((n, n), dtype=float)
        sip[j, i] = v
        return sip

    def load_flatfield(self, **kwargs):
        """
        Method to load a flat field

        Parameters
        ----------
        kwargs : dict, optional
            Optional variables passed to `load_flatfield`

        Returns
        -------
        flat : `np.ndarray`
            the 2d flat field.

        """
        flat = self.config.load_flatfield(**kwargs)

        return flat

    def items(self):
        """
        An iterator to emulate the dict-like behavior over the spectral orders

        Returns
        -------
        2-tuple : ordername and `Order` class
        """
        yield from self.config.items()

    def keys(self):
        """
        Method to emulate the dict-like behavior over the spectral orders

        Returns
        -------
        list : the order names
        """
        return self.config.keys()

    def values(self):
        """
        Method to emulate the dict-like behavior over the spectral orders

        Returns
        -------
        list : the `Order` objects
        """
        return self.config.values()

    def __getitem__(self, k):
        """
        Method to emulate the dict-like behavior over the spectral orders
        """
        return self.config[k]

    def __iter__(self):
        """
        An iterator over the spectral orders
        """
        yield from self.config

    def drizzle(self, x0, y0, order, wav, **kwargs):
        """
        Clip a closed polygon against a pixel grid and collect fractional
        pixel areas in the new grid.

        Parameters
        ----------
        x0 : list, tuple, `np.ndarray`
            A list of x-coordinate polygon indices

        y0 : list, tuple, `np.ndarray`
            A list of y-coordinate polygon indices

        order : str
            the name of the spectral order

        wav : list, tuple, `np.ndarray`
            A list of wavelengths

        Returns
        -------
        x : `np.ndarray` (dtype==int)
            The x pixel coordinates in the new grid.  Will be an integer

        y : `np.ndarray` (dtype==int)
            The y pixel coordinates in the new grid.  Will be an integer

        lam : `np.ndarray` (dtype==uint16)
            The wavelength indices

        area : `np.ndarray` (dtype == float)
            The fractional pixel areas

        """

        # check that pixel is disperable
        if not self.config.pom(np.average(x0), np.average(y0)):
            return [], [], [], []

        # disperse the pixel
        xg, yg = self.config.disperse(x0, y0, order, wav)

        # truncate pixels to be in the grid (move this into pypolyclip?)
        xg = np.clip(xg - 0.5, 0, self.naxis[0] - 1)
        yg = np.clip(yg - 0.5, 0, self.naxis[1] - 1)

        # clip against the pixel grid
        x, y, area, slices = clip_multi(xg, yg, self.naxis)

        # make wavelength indices
        lam = np.empty_like(x, dtype=np.uint16)
        for i, s in enumerate(slices):
            lam[s] = i

        return x, y, lam, area


class InstrumentConfig(dict):
    """
    Class to configure an instrument

    inherits from dict.  The key, value pairs are for the detector
    name and the `DetectorConfig` object
    """

    def __init__(self, telescope, instrument, disperser, subarray=False, **kwargs):
        """
        Initializer

        Parameters
        ----------
        telescope : str
            The name of the telescope

        instrument : str
            The name of the instrument

        disperser : str or 2-tuple
            A variable defining the type of disperser:

            If a `str`, then the blocking is assumed as `None`
            If a 2-tuple, then assumed they are the disperser and blocking name

            This is used to create a `Disperser` object.

        kwargs : dict, optional
            Keywords to configure `DetectorConfig` objects


        Notes
        -----
        Read a yaml file from the `slitlessutils_config` directory

        """
        self.filename = os.path.join(Config().refpath, 'instruments',
                                     f'{telescope}_{instrument}.yaml'.lower())

        # load the yaml config file
        with open(self.filename) as f:
            data = yaml.safe_load(f)

        # extract relevant data from the yaml dict
        self.telescope = data['telescope']
        self.instrument = data['instrument']
        self.bunit = data['bunit']
        self.suffix = data['suffix']
        self.subarray = subarray

        self.path = data['path']
        self.header = data.get('header', {})
        self.siaf = SIAF(*data['siaf'].values())

        if isinstance(disperser, str):
            blocking = None
        elif isinstance(disperser, (tuple, list)) and len(disperser) == 2:
            disperser, blocking = disperser
        else:
            raise TypeError(f"disperser is invalid: {disperser}")

        # make some more stable error checking
        try:
            dispdata = data['dispersers'][disperser]
        except KeyError:
            msg = f"disperser {disperser} is not found in reference manifest"
            LOGGER.error(msg)
            raise KeyError(msg)

        try:
            d = dispdata[blocking]
        except KeyError:
            msg = f"blocking {blocking} is not found in reference manifest"
            LOGGER.error(msg)
            raise KeyError(msg)

        self.disperser = load_disperser(d['extraction'], name=disperser, blocking=blocking)
        self.tabulator = load_disperser(d['tabulation'], name=disperser, blocking=blocking)

        # self.disperser = load_disperser(disperser, blocking, **d)
        self.backgrounds = d.get('background', {})
        for detname, detdata in data['detectors'].items():
            self[detname] = DetectorConfig(detname, detdata, self.disperser,
                                           data['path'], **kwargs)

    @classmethod
    def from_fitsfile(cls, filename, **kwargs):
        """
        Classmethod to load an `InstrmentConfig` from a fits file (ie. the
        content in the primary header)

        Parameters
        ----------
        filename : str
            Full path to the fits file

        kwargs : dict, optional
            Variables that are passed to the initializer

        Returns
        -------
        inscnf : `InstrumentConfig`
            The instrument configuration object
        """
        with fits.open(filename, mode='readonly') as hdul:
            h0 = hdul[0].header

            # check that the OBSTYPE is a spectroscopic image
            if h0['OBSTYPE'] == 'SPECTROSCOPIC':
                tel = h0['TELESCOP']
                ins = h0['INSTRUME']
                subarray = h0.get('SUBARRAY', False)

                # do something special when loading each type of instrument
                if ins == 'WFC3':
                    ins += h0['DETECTOR']
                    disperser = h0['FILTER']

                elif ins == 'ACS':
                    ins += h0['DETECTOR']
                    disperser = h0['FILTER1']

                elif ins == 'NIRISS':
                    disperser = (h0['FILTER'], h0['PUPIL'])
                    kwargs['fwcpos'] = h0['FWCPOS']
                else:
                    raise NotImplementedError(f"Unsupported: {tel = } {ins = }")

                # instantiate the object
                insconf = cls(tel, ins, disperser, subarray=subarray, **kwargs)

                # update some things because subarrays
                detnames = list(insconf.keys())
                for detname in detnames:
                    ext = insconf[detname].extensions['science'].extension

                    # gotta check this to deal with subarrays with instruments that
                    # have multiple detectors (like ACS/WFC or WFC3/UVIS)
                    if ext in hdul:
                        h = fits.getheader(filename, ext)

                        if subarray:
                            # if this is a subarray-ed image, then update dimensions
                            if insconf[detname].naxis[0] != h['NAXIS1']:
                                insconf[detname].naxis[0] = h['NAXIS1']
                            if insconf[detname].naxis[1] != h['NAXIS2']:
                                insconf[detname].naxis[1] = h['NAXIS2']
                            if insconf[detname].crpix[0] != h['CRPIX1']:
                                insconf[detname].crpix[0] = h['CRPIX1']
                            if insconf[detname].crpix[1] != h['CRPIX2']:
                                insconf[detname].crpix[1] = h['CRPIX2']

                    else:
                        # remove extensions for a subarray
                        del insconf[detname]

            else:
                insconf = None
        return insconf

    def background_filename(self, key):
        """
        A method to return the filename of the background image(s)

        Parameters
        ----------
        key : str
            The type of background to retrieve

        Returns
        -------
        filename : str
            The full path to the requested sky image.  If key does not exist,
            then returns `None`
        """
        back = self.backgrounds.get(key)
        if back:
            back = os.path.join(Config().refpath, 'instruments',
                                self.path, back)
        return back

    def npixels(self):
        """
        Method to get the total number of pixels in a single file

        Returns
        -------
        npix : int
            The total number of pixels
        """
        npix = sum(np.prod(v.naxis) for v in self.values())
        return npix

    def equivalent(self, telescope, instrument, disperser):
        """
        Method to test equivalence, which is defined here as if this
        instance of an `InstrumentConfig` is for the same telescope,
        instrument, and disperser/blocking

        The confounding issue is that the JWST/NIRISS and NIRCam grisms
        will be the same disperser (e.g. GR150 for NIRISS), but will be in
        different orientation  (e.g. R vs C).  Nevertheless they are (may be)
        equivalent for spectral extraction depending on the purpose.

        Parameters
        ----------
        telescope : str
            The name of the telescope

        instrument : str
            The name of the instrument

        disperser : 2-tuple, str, or `Disperser`
            If a 2-tuple, then the elements are the disperser and blocking
            filter's name, respectively.  If a str, then only the disperser
            name is checked.  If a `Disperser` object, then all attributes
            are checked.  In the first two cases, the disperser is checked
            against the "trimmed" name, which removes the last character.

        Returns
        -------
        logic : bool
            A flag if they're equivalent
        """

        if isinstance(disperser, tuple) and len(disperser) == 2:
            return (self.telescope == telescope) and \
                (self.instrument == instrument) and \
                (self.disperser.trimmed_name == disperser[0]) and \
                (self.disperser.blocking == disperser[1])

        elif isinstance(disperser, str):
            return (self.telescope == telescope) and \
                (self.instrument == instrument) and \
                (self.disperser.trimmed_name == disperser)

        elif isinstance(disperser, Disperser):
            return (self.telescope == telescope) and \
                (self.instrument == instrument) and \
                (self.disperser.trimmed_name == disperser.trimmed_name) and \
                (self.blocking == disperser.blocking)
        else:
            return False

    def make_header(self, targname=None):
        """
        Method to create a primary fits header from scratch

        Parameters
        ----------
        targname : str or `None`
            The name of the target.  If `None` then will be interpreted
            as blank.  Default is `None`.

        Returns
        -------
        phdr : `astropy.io.fits.Header`
            The fits header

        """

        phdr = fits.Header()

        # for k,v in self.header.items():
        #    phdr[k]=tuple(v)

        phdr['TELESCOP'] = (self.telescope, 'telescope used to acquire the data')

        # process each instrument separately
        if self.instrument in ("WFC3IR", 'WFC3UVIS'):
            phdr['INSTRUME'] = (self.instrument[:4],
                                'identifier for instrument used to acquire data')
            phdr['DETECTOR'] = (self.instrument[4:], 'detector in use: UVIS or IR')
            phdr['FILTER'] = (self.disperser.name, 'element selected from filter wheel')
        elif self.instrument in ('ACSWFC', 'ACSSBC'):
            phdr['INSTRUME'] = (self.instrument[:3],
                                'identifier for instrument used to acquire data')
            phdr['DETECTOR'] = (self.instrument[3:], 'detector in use: WFC, HRC, or SBC')
            phdr['FILTER1'] = (self.disperser.name, 'element selected from filter wheel 1')
            phdr['FILTER2'] = ('CLEAR2L', 'element selected from filter wheel 2')
        elif self.instrument == 'NIRISS':
            phdr['INSTRUME'] = ('NIRISS', 'Instrument used to acquire the data')
            phdr['DETECTOR'] = ('NIS', 'Name of detector used to acquire the data')
            phdr['FILTER'] = (self.disperser.name, 'element selected from filter wheel')
            phdr['PUPIL'] = (self.disperser.blocking, 'element selected from pupil wheel')
        else:
            raise NotImplementedError(f"{self.telescope}_{self.instrument}")

        phdr['OBSTYPE'] = ('SPECTROSCOPIC', 'observation type - imaging or spectroscopic')
        headers.add_stanza(phdr, 'Instrument Configuration Information',
                           before='TELESCOP')

        if not isinstance(targname, str):
            targname = ' ' * 8
        phdr['TARGNAME'] = (targname, "proposer's target name")
        headers.add_stanza(phdr, 'Target Information', before='TARGNAME')

        # put some timing information in
        dt = datetime.now()

        exptime = self[list(self.keys())[0]].noise.time
        delta = dt - MJDREF
        mjd0 = delta.days + (delta.seconds + delta.microseconds / 1e6) / 86400.
        mjd1 = mjd0 + exptime / 86400.

        phdr['DATE-OBS'] = (dt.strftime('%Y-%m-%d'), 'UT date of start of observation (yyyy-mm-dd)')
        phdr['TIME-OBS'] = (dt.strftime('%H:%M:%S'), 'UT time of start of observation (hh:mm:ss)')
        phdr['EXPSTART'] = (mjd0, 'exposure start time (Modified Julian Date)')
        phdr['EXPEND'] = (mjd1, 'exposure end time (Modified Julian Date)')
        phdr['EXPTIME'] = (exptime, 'exposure duration (seconds)')
        phdr['EXPFLAG'] = ('NORMAL', 'Exposure interruption indicator')

        headers.add_stanza(phdr, 'Exposure Information', before='DATE-OBS')

        self.disperser.update_header(phdr)
        self.siaf.update_header(phdr)

        return phdr

    def make_filename(self, dataset):
        """
        Method to create a filename

        Parameters
        ----------
        dataset : str
             The name of the dataset

        Returns
        -------
        filename : str
             The name of the fits file with a suffix added

        """
        filename = f'{dataset}_{self.suffix}.fits'
        return filename
