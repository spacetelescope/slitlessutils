import os
import warnings

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS, Sip
from astropy.wcs import utils as wcsutils

from ....logger import LOGGER
from ..config import InstrumentConfig


class WFSSDetector:
    """
    Class to contain information for a single WFSS image detector.

    """

    def __init__(self, filename, detconf, bunit='e-/s'):
        """
        Initializer method

        Parameters
        ----------
        filename : str
            full path to a valid fits file that contains this data

        detconf :
            Configuration object for this detector

        bunit : str, optional
            the BUNIT of the data.  Default is 'e-/s'

        """

        # save the filename
        self.filename = filename

        # create an empty WCS to get started
        self.wcs = WCS(relax=True)
        self.wcs._naxis = list(detconf.naxis)  # NB: MUST BE A LIST!

        # do we add SIP?
        if detconf.sip:
            self.wcs.wcs.ctype = ['RA---TAN-SIP', 'DEC--TAN-SIP']
            self.wcs.sip = Sip(*detconf.sip)
        else:
            self.wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']

        # record the CRPIX & CDELT
        self.wcs.wcs.crpix = detconf.crpix
        self.wcs.wcs.cdelt = detconf.scale / 3600.

        # save the config object
        self.config = detconf
        self.bunit = bunit

    def get_pa(self, degrees=True, limit=2.0):
        """
        Method to compute the position angle, defined as the angle between
        +y and N axes, at the CRPIX from a WCS object

        Parameters
        ----------
        degrees : bool, optional
            Flag to specify the units.  Default is True

        limit : float, optional
            Limit between two possible angles.  If large, then this is
            likely due to skew.  Default is 2 deg.

        Returns
        -------
        theta : float
            The angle
        """

        # get the handedness of the coordinates
        sgn = np.sign(np.linalg.det(self.wcs.wcs.piximg_matrix))

        # compute the two angles (+y to N and -x to E)
        angx = np.arctan2(sgn * self.wcs.wcs.piximg_matrix[0, 1],
                          sgn * self.wcs.wcs.piximg_matrix[0, 0])
        angy = np.arctan2(-self.wcs.wcs.piximg_matrix[1, 0],
                          +self.wcs.wcs.piximg_matrix[1, 1])

        # check the angular difference, mindful of the wrapping at on [0,2pi)
        dang = (angx - angy + np.pi) % (2 * np.pi) - np.pi

        # issue a quick warning
        if np.abs(dang) < np.radians(limit):
            msg = f'X and Y axes rotations differ by more than {limit} deg.'
            LOGGER.warning(msg)

        if degrees:
            return np.degrees(angy)
        else:
            return angy

    def get_pixscl(self):
        """
        Method to compute the pixel scale (in arcsec/pix) at the CRPIX
        from a WCS object

        Returns
        -------
        px : float
            pixel scale in x in arcsec/pix

        py : float
            pixel scale in y in arcsec/pix

        Notes
        -----
        Often the x-pixel scale will be negative --- this reflects the
        N-up, E-left in usual astronomical coordinate systems.
        """

        # get the handedness of the coordinates
        sgn = np.sign(np.linalg.det(self.wcs.wcs.piximg_matrix))

        # a dummy variable
        c2 = self.wcs.wcs.piximg_matrix**2

        # the pixel scales
        px = sgn * np.sqrt(np.sum(c2[:, 0])) * 3600.
        py = np.sqrt(np.sum(c2[:, 1])) * 3600.

        return px, py

    def __str__(self):
        return self.name

    @property
    def extensions(self):
        return self.config.extensions

    @property
    def orders(self):
        return self.config.config.keys()

    @property
    def pixelarea(self):
        return wcsutils.proj_plane_pixel_area(self.wcs)

    @property
    def naxis(self):
        return self.wcs._naxis

    @property
    def shape(self):
        return self.wcs._naxis[::-1]

    @property
    def name(self):
        return self.config.name

    def npixels(self):
        """
        Method to return the number of pixels in this detector

        Returns
        -------
        npix : str
           number of pixels
        """
        npix = self.wcs._naxis[0] * self.wcs._naxis[1]
        return npix

    def relative_pixelarea(self, x, y):
        """
        Method to compute the pixel area relative to the pixel area
        at the CRPIX.  This is equivalent to the pixel-area map (PAM)

        Computes the Jacobian of the distortion solution on the fly, therefore
        if no distortion is present, then returns all unity.


        Parameters
        ----------
        x : float, int, or `np.ndarray`
            the x-coordinates in the detector frame to compute

        y : float, int, or `np.ndarray`
            the y-coordinates in the detector frame to compute

        Returns
        -------
        area : float or `np.ndarray`
            the relative pixel area --- if an `np.ndarray`, then the
            dtype will be `float`

        Notes
        -----
        The input coordinates (x,y) must have the same shape, but their
        dtype does not matter.

        """

        if self.wcs.sip:
            # if SIP is present, then compute the Jacobian
            dx = x - (self.wcs.wcs.crpix[0] - 1)
            dy = y - (self.wcs.wcs.crpix[1] - 1)

            # Uggh...computing these derivatives is ugly, because
            # things are stored as 2d matrices.

            # compute the derivatives w.r.t. the A
            dadx = np.ones_like(dx, dtype=float)
            dady = np.zeros_like(dx, dtype=float)
            for i, j in zip(*np.where(self.wcs.sip.a != 0)):
                if i > 0:
                    dadx += self.wcs.sip.a[i, j] * i * dx**(i - 1) * dy**j
                if j > 0:
                    dady += self.wcs.sip.a[i, j] * j * dx**i * dy**(j - 1)

            # compute the derivatives w.r.t. the B
            dbdx = np.zeros_like(dx, dtype=float)
            dbdy = np.ones_like(dx, dtype=float)
            for i, j in zip(*np.where(self.wcs.sip.b != 0)):
                if i > 0:
                    dbdx += self.wcs.sip.b[i, j] * i * dx**(i - 1) * dy**j
                if j > 0:
                    dbdy += self.wcs.sip.b[i, j] * j * dx**i * dy**(j - 1)

            # compute the jacobian:
            area = np.abs(dadx * dbdy - dady * dbdx)
        else:
            # if no distorion, then return all unity
            area = np.ones_like(x, dtype=float)

        return area

    def set_pointing(self, ra, dec, orientat, refsiaf):
        """
        Method to set the pointing to the WCS structure.

        This is primarily only useful when setting up a simulated WFSS image

        Parameters
        ----------
        ra : float or int
            The RA for this detector (will become the CRVAL1)

        dec : float or int
            The Dec for this detector (will become the CRVAL2)

        orientat : float or int
            The position angle of the detector in deg (is angle East of North)

        refsiaf : `su.core.wfss.config.SIAF`
            The SIAF data that describes the overall instrumental pointing.
            Used to find the offset from the instrument pointing to this
            detector's pointing.

        Returns
        -------
        None
        """

        # orientat = PA_APER for JWST
        pa = np.radians(orientat)
        cs = np.cos(pa)
        sn = np.sin(pa)

        self.wcs.wcs.pc = [[-cs, sn], [sn, cs]]
        self.wcs.wcs.crval = self.config.siaf.compute_crvals(ra, dec, orientat,
                                                             refsiaf)

    def set_wcs(self, wcs):
        """
        A method to set the WCS

        Parameters
        ----------
        wcs : `astropy.wcs.WCS`
            the WCS to set to this detector
        """
        self.wcs = wcs.copy()

    def make_header(self, imtype):
        """
        Method to make a fits header

        Parameters
        ----------
        imtype : str
            The type of the image: 'science', 'uncertainty', or 'dataquality'

        Returns
        -------
        h : `astropy.io.fits.Header`
            The fits header

        """

        if imtype in self.extensions:

            h = self.wcs.to_header(relax=True)
            for k, v in self.config.header.items():
                h[k] = tuple(v)
            self.extensions[imtype].update_header(h)
            self.config.config.pom.update_header(h)
        else:
            h = fits.Header()

        return h

    def make_HDUs(self, sci, addnoise=True):
        """
        Method to create header-data units (HDUs) for the science,
        uncertainty, and data-quality.

        Parameters
        ----------
        sci : `np.ndarray`
            Science image --- a two-dimensional array

        addnoise : bool, optional
            Flag to add noise to the image. Default is True.

        Returns
        -------
        scihdu : `astropy.io.fits.ImageHDU`
            The HDU containing the science image and header.  If
            addnoise==False, then this is effectively the same as the
            input with units changed based on BUNIT.

        unchdu : `astropy.io.fits.ImageHDU`
            The HDU containing the uncertainty image and header

        dqahdu : `astropy.io.fits.ImageHDU`
            The HDU containing the data-quality image and header

        """

        hsci = self.make_header('science')
        hunc = self.make_header('uncertainty')
        hdqa = self.make_header('dataquality')

        if addnoise:
            sci, unc = self.config.noise(sci)
            self.config.noise.update_header(hsci)

        else:
            unc = np.ones_like(sci)

        # change the units to match BUNIT
        hsci['BUNIT'] = (self.bunit, 'brightness unit')
        hunc['BUNIT'] = (self.bunit, 'brightness unit')
        hdqa['BUNIT'] = ('UNITLESS', 'brightness unit')
        if self.bunit.lower() in ('electron', 'electrons', 'e', 'e-'):
            sci *= self.config.noise.time
            unc *= self.config.noise.time

        scihdu = fits.ImageHDU(data=self.extensions['science'].retype(sci),
                               header=hsci)
        unchdu = fits.ImageHDU(data=self.extensions['uncertainty'].retype(unc),
                               header=hunc)
        dqahdu = fits.ImageHDU(
            data=np.zeros_like(
                sci,
                dtype=self.extensions['dataquality'].dtype),
            header=hdqa)

        return scihdu, unchdu, dqahdu

    def readfits(self, imtype, header=False):
        """
        Method to read and return data from a fits file

        Parameters
        ----------
        imtype : str
             The type of the image to read: 'science', 'uncertainty' or
             'dataquality'

        header : bool, optional
             Flag to additionally read the header.  Default is False

        Returns
        -------
        img : `np.ndarray`
             The image (a two-dimensional array)

        hdr : `astropy.io.fits.Header`, optional
             The fits header.  Default is to not return the header

        """

        if imtype in self.extensions:
            return self.extensions[imtype].readfits(self.filename,
                                                    header=header)

            # args=self.extensions[imtype].readfits(self.filename,header=header)
            #
            # if header:
            #    if self.bunit in ('electrons','e','e-','electron'):
            #        return args[0]/self.config.noise.time,args[1]
            #    else:
            #        return args
            # else:
            #    if self.bunit in ('electrons','e','e-','electron'):
            #        return args/self.config.noise.time,args[1]
            #    else:
            #        return args
        elif imtype == 'flatfield':
            return fits.getdata(self.flatfield, header=header)
        else:
            pass

    def primaryheader(self):
        """
        Method to read the primary header

        Returns
        -------
        phdr : `astropy.io.fits.Header`
            The primary header.
        """
        phdr = fits.getheader(self.filename, 0)
        return phdr

    def headfits(self, imtype):
        """
        Method to read a fits header from a file

        Parameters
        ----------
        imtype : str
            The type of the image: 'science', 'uncertainty', or 'dataquality'

        Returns
        -------
        h : `astropy.io.fits.Header`
            The fits header
        """
        if imtype in self.extensions:
            h = self.extensions[imtype].headfits(self.filename)
        else:
            h = fits.Header()
        return h

    def xy2xy(self, x, y, wcs, forward=True):
        """
        Method to transform (x,y) pairs from one WCS frame to a different
        WCS frame

        Parameters
        ----------
        x : int, float, `np.ndarray`
           The x-coordinates

        y : int, float, `np.ndarray`
           The y-coordinates

        wcs : astropy.wcs.WCS
           The other WCS to transform to

        forward : bool, optional
           A flag controlling the direction of the transform.

           -> If True, then transform *FROM* this image *TO* this WCS.
           -> If False, then transform *TO* this image *FROM* this WCS

        Returns
        -------
        X : float, `np.ndarray`
           the transformed x-coordinates

        Y : float, `np.ndarray`
           the transformed y-coordinates

        Notes
        -----
        The input `x` and `y` variables must have the same shape, and then
        the output `X` and `Y` will have that shape
        """

        if forward:
            #X, Y = wcs.all_world2pix(*self.wcs.all_pix2world(x, y, 0), 0)
            X, Y = wcs.wcs_world2pix(*self.wcs.wcs_pix2world(x, y, 0), 0)
        else:
            #X, Y = self.wcs.all_world2pix(*wcs.all_pix2world(x, y, 0), 0)
            X, Y = self.wcs.wcs_world2pix(*wcs.wcs_pix2world(x, y, 0), 0)
        return X, Y

    def ad2xy(self, a, d):
        """
        Method to transform (RA,Dec) pairs to (x,y) pairs.

        Parameters
        ----------
        a : int, float, `np.ndarray`
           The RA coordinates

        d : int, float, `np.ndarray`
           The Dec coordinates

        Returns
        -------
        x : float, `np.ndarray`
           The x-coordinates

        y : float, `np.ndarray`
           The y-coordinates


        Notes
        -----
        The input `a` and `d` variables must have the same shape, and the
        output `x` and `y` will have that shape.
        """
        x, y = self.wcs.all_world2pix(a, d, 0)
        return x, y

    def xy2ad(self, x, y):
        """
        Method to transform (x,y) pairs to (RA,Dec) pairs.

        Parameters
        ----------
        x : int, float, `np.ndarray`
           The x-coordinates

        y : int, float, `np.ndarray`
           The y-coordinates

        Returns
        -------
        a : float, `np.ndarray`
           The RA coordinates

        d : float, `np.ndarray`
           The Dec coordinates

        Notes
        -----
        The input `x` and `y` variables must have the same shape, and the
        output `a` and `d` will have that shape.
        """
        a, d = self.wcs.all_pix2world(x, y, 0)
        return a, d


class WFSS(dict):
    """
    Class to contain all the information for a single WFSS file

    inherits from dict.  The keyword,value pairs are for each detector will
    be detectorname,`WFSSDetector` respectively.
    """

    def __init__(self, filename, filetype, insconf):
        """
        Initializer

        Parameters
        ----------
        filename : str
            The name of the input file (full path)

        filetype : str
            The type of the file ('observed' or 'simulated')

        insconf : `slitlessutils.core.wfss.config.InstrumentConfig`
            The instrument configuration object

        """

        self.filename = filename
        self.filetype = filetype
        self.config = insconf

        for detname, detconf in self.config.items():
            self[detname] = WFSSDetector(self.filename, detconf,
                                         bunit=insconf.bunit)

        self.visit = ''
        self.orbit = 1
        self.subarray = False

    def __str__(self):
        s = [f' WFSS file: {self.filename}',
             # f' \033[39;4;1mWFSS file: {self.filename}\033[00m',
             f' file type: {self.filetype}',
             f' telescope: {self.telescope}',
             f'instrument: {self.instrument}',
             f' disperser: {self.disperser.name}',
             f'  blocking: {self.disperser.blocking}']

        return '\n'.join(s)

    def extensions(self):
        """
        An iterator to loop over the detectors and all possible extensions

        Returns
        -------
        detname : str
            The name of the detector

        extname : str
            The name of the extension

        extdata : `np.ndarray`
            The data from this extension/detector
        """
        for detname, detdata in self.items():
            for extname, extdata in detdata.config.extensions.items():
                yield (detname, extname), extdata

    def get_pa(self, **kwargs):
        """
        Method to get the average PA over all the detectors

        Parameters
        ----------
        kwargs : dict, optional
            See `WFSSDetector.get_pa()` for the options

        Returns
        -------
        pa : float
            The average PA
        """
        n = len(self)
        pas = np.zeros(n, dtype=float)
        for i, det in enumerate(self.values()):
            pas[i] = det.get_pa(**kwargs)

        return np.average(pas)

    def get_pixscl(self):
        """
        Method to get the average pixel scales over all the detectors

        Returns
        -------
        px : float
            The average x-axis pixel scale

        py : float
            The average y-axis pixel scale
        """

        n = len(self)
        pxs = np.zeros(n, dtype=float)
        pys = np.zeros(n, dtype=float)
        for i, det in enumerate(self.values()):
            pxs[i], pys[i] = det.get_pixscl()

        return np.average(pxs), np.average(pys)

    def get_parameters(self):
        """
        Method to read the dispering information

        Returns
        -------
        disperser : `su.core.wfss.config.Disperser`
            the information on the dispersive element
        """

        conf = InstrumentConfig.from_fitsfile(self.filename)
        return conf.disperser

    @property
    def dataset(self):
        filename = os.path.basename(self.filename)
        tokens = filename.split(f'_{self.config.suffix}')
        return tokens[0]

    @property
    def telescope(self):
        return self.config.telescope

    @property
    def instrument(self):
        return self.config.instrument

    @property
    def disperser(self):
        return self.config.disperser

    @property
    def blocking(self):
        return self.config.disperser.blocking

    @property
    def siaf(self):
        return self.config.siaf

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
        return self.config.background_filename(key)

    @classmethod
    def simulated(cls, telescope, instrument, dataset, ra, dec, orientat, dispname,
                  **kwargs):
        """
        Classmethod to load simulated WFSS data

        Parameters
        ----------
        telescope : str
            The name of the telescope

        instrument : str
            The name of the instrument

        dataset : str
            The name of the dataset

        ra : float
            The center of the field in RA in deg

        dec : float
            The center of the field in Dec in deg

        orientat : float
            The rotation of the field in deg.  Angle is measured East of North

        dispname : str
            The name of the dispersive element

        kwargs : dict, optional
             Optional keywords that are passed to
             `slitlessutils.core.wfss.config.InstrumentConfig`.  A key
             variable one might add is `blocking=<str>` to enable
             blocking filters (a la JWST)

        Returns
        -------
        obj : `WFSS`
             The WFSS data structure

        """
        insconf = InstrumentConfig(telescope, instrument, dispname, **kwargs)
        obj = cls(insconf.make_filename(dataset), 'simulated', insconf)

        # compute somethings
        obj.pa_v3 = insconf.siaf.compute_pav3(orientat)    # PA_V3

        for detname, detdata in obj.items():
            detdata.set_pointing(ra, dec, orientat, insconf.siaf)

        return obj

    @classmethod
    def observed(cls, filename, **kwargs):
        """
        Classmethod to load observed WFSS data

        Parameters
        ----------
        filename : str
            The fullpath to the fits file

        kwargs : dict, optional
            keywords passed to
            `slitlessutils.core.wfss.config.InstrumentConfig`.

        Returns
        -------
        obj : `WFSS`
            The WFSS object
        """
        insconf = InstrumentConfig.from_fitsfile(filename, **kwargs)

        obj = cls(filename, 'observed', insconf)

        with fits.open(filename, mode='readonly') as hdul:
            for det in obj.values():
                with warnings.catch_warnings():
                    warnings.filterwarnings("ignore", 'FITSFixedWarning')
                    wcs = WCS(det.headfits('science'), hdul, relax=True)

                det.set_wcs(wcs)

        # parse the zeroth header to get some global data
        phdr = fits.getheader(obj.filename, ext=0)
        tel = phdr.get('TELESCOP')
        if tel == 'HST':

            # test for SUBARRAY
            obj.subarray = phdr.get('SUBARRAY', False)
            if obj.subarray:
                LOGGER.knownissue(f'Slitlessutils does not fully support subarray data: {filename}')

            # get the VISIT info
            line = phdr.get('LINENUM')
            if line:
                s = line.split('.')
                obj.visit = s[0]
                obj.orbit = int(s[1])
        elif tel == 'JWST':
            obj.visit = phdr['VISIT']
            obj.orbit = 0
            obj.subarray = False
        else:
            LOGGER.error(f'Telescope ({tel}) is not found to get visit')

        return obj
