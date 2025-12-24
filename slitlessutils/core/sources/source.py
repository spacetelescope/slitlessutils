import os

import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.modeling import fitting, models
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS
from astropy.wcs import utils as wcsutils
from skimage.segmentation import expand_labels

from ...logger import LOGGER
from ..utilities import headers, indices
from .dispersedregion import DispersedRegion


class Source(list):
    """
    Class to hold a single spectral source, which in turn may contain any
    number of spectral regions that have a unique spectrum.

    inherits from list.  Each element of the list is a `DispersedRegion`,
    which controls the unique spectrum.
    """

    # maximum number of spectral regions
    MAXSPECREG = 1000

    def __init__(self, segid, img, seg, hdr, reg=None, zeropoint=26., grpid=0,
                 local_back=True, backsize=5, nsig=(5, 5),
                 whttype='pixels', profile='gaussian', negfunc='positivity',
                 epsilon=1e-9):
        """
        Initializer

        Parameters
        ----------
        img : `np.ndarray`
            A two-dimensional array from the direct image

        seg : `np.ndarray`
            A two-dimensional array from the segmentation image

        reg : `np.ndarray`, optional
            If an np.array, then this is a "segmentation map" like object,
            then these pixels get aggregated.  If None, then creates a flat
            region image.  Default is None

        hdr : `astropy.io.fits.Header`
            A fits header for the passed direct image

        local_back : bool, optional
            Flag to subtract the local background in the direct image.  In
            general, the direct image pixel values are factors in the
            extraction weights, and so if there is an additive pedestal
            constant, then the weights will be very suspect.  Therefore,
            one should try to remove the local background in the direct
            image.  Default is True

        backsize : int, optional
            The size (in pixels) of the local background annulus.
            Default is 5

        nsig : 2-tuple, int, or float, optional
            The number of sigma to clip the sky background.  If the value
            is a tuple, then the numbers are interpreted as (low,high)
            values. If it is a single value, then will be both low and high
            clipping values.  Default is (5.,5.).

        whttype : str, optional
            Way to control how the direct image values are processed into
            the weights.  If equal to:

            'pixels'       Use the absolute value of the pixels
            'fitprofile'   Fit a profile with some analytic profile, as
                           specified by the `profile` variable.

        profile : str, optional
            Name of the profile to fit.  See `astropy.modeling.models`.
            This is only used if ```whttype='fitprofile'```.  Presently,

            gaussian:  will fit a 2d-Gaussian
            sersic:  will fit a 2d-Sersic

            Default is 'gaussian'
              NOTE: this functionality is not intended to provide reliable
                    estimates for the profile's parameters, but only a
                    smooth description of the source's spatial profile. As
                    such, the parameters are not available for output.

        zeropoint : float or int, optional
            The photometric zeropoint in the AB mag system for the direct
            image. This is used in two ways:
            1) to compute the magnitude through the segmentation region
            2) to normalize an input SED (such as when simulating)
            Default is 26.

        grpid : int, optional
            The group identifier.  A 'group' is a collection of sources
            that will have spectral overlap in this series of WFSS images.
            see `su.core.modules.Group` for more information. Default is 0.

        epsilon : float, optional
            The critical value for forcing positivity.  Any value below this
            will be re-mapped onto the positive domain.  Default is 1e-9.
            See the IDL implementation for more info.
            https://idlastro.gsfc.nasa.gov/ftp/pro/image/positivity.pro

        Notes
        ----
        This is a cornerstone object to slitlessutils.

        """

        self.wcs = WCS(hdr, relax=True)
        self.ltv = [hdr.get('LTV1', 0.), hdr.get('LTV2', 0.)]

        # record some things
        self.segid = segid
        self.grpid = grpid
        # self.whttype = whttype
        self.backsize = max(backsize, 0)
        self.local_back = local_back

        # set the spectral parameters
        self.specpars = {'lamb0': None, 'lamb1': None, 'dlamb': None}

        # test the number of sigma for the sky
        if isinstance(nsig, (float, int)):
            nsig = (nsig, nsig)
        self.nsig = nsig

        # find the pixels for this segmentation
        y, x = np.where(seg == self.segid)

        # check for a valid image
        self.npixels = len(x)
        if self.npixels > 0:

            # crop the image
            # get a bounding box
            x0 = max(np.amin(x) - self.backsize, 0)
            y0 = max(np.amin(y) - self.backsize, 0)
            x1 = min(np.amax(x) + self.backsize, img.shape[1] - 1)
            y1 = min(np.amax(y) + self.backsize, img.shape[0] - 1)

            # cut out regions
            subimg = img[y0:y1, x0:x1]
            subseg = seg[y0:y1, x0:x1]
            # subwht = wht[y0:y1, x0:x1]

            # compute and subtract the local sky background:
            if self.local_back:

                # grow the region to find the local sky background
                # but use a mask in case it results in empty array
                subseg2 = expand_labels(subseg, distance=self.backsize)
                submsk = ~np.logical_or(subseg2 != self.segid, subseg != 0)
                try:
                    ave, med, sig = sigma_clipped_stats(subimg, mask=submsk,
                                                        sigma_lower=self.nsig[0],
                                                        sigma_upper=self.nsig[1])
                except BaseException:
                    ave = 0.0

                self.background = ave
                subimg -= self.background
            else:
                self.background = 0.0

            # convert the image into weights
            w = self.compute_weights(subimg, x - x0, y - y0, whttype=whttype, profile=profile,
                                     negfunc=negfunc, epsilon=1e-9)

            # remove pixels with zero weights
            g = np.where(w > 0)[0]
            self.npixels = len(g)
            if self.npixels == 0:
                return
            x = x[g]
            y = y[g]
            w = w[g]

            # compute the centeroids
            xyc = (np.average(x, weights=w), np.average(y, weights=w))
            self.xyc = self.image_coordinates(*xyc)
            adc = self.wcs.all_pix2world(*xyc, 0)
            self.adc = (float(adc[0]), float(adc[1]))

            # compute the area of the object
            self.area = self.npixels * self.pixelarea

            #
            # compute some fluxes:
            # flux = instrumental flux (in e-)
            # fnu  = physical flux in (erg/s/cm2/Hz)
            # mag  = AB mag.
            #
            # NOTA BENE: norm and flux would be exactly the same iff
            #            whttype=='pixels', else they'll be similar
            #
            self.flux = np.sum(img[y, x])
            self.fnu = self.flux * 10.**(-0.4 * (zeropoint + 48.6))
            if self.flux < 0:
                self.mag = np.nan
            else:
                self.mag = -2.5 * np.log10(self.flux) + zeropoint

            # put coordinates back on original footprint
            # x,y=self.image_coordinates(x,y,dtype=np.int)

            # Now parse the source for the DispersedRegions.
            if isinstance(reg, np.ndarray) and reg.shape == seg.shape:
                # have a valid region image, so use it
                r = reg[y, x].astype(int)
            elif self.segid < 0:
                # a negative SEGID means checkerboard REGION
                r = np.arange(1, len(x) + 1, dtype=int)
            else:
                # if everything is invalid, then set to a single value
                r = np.ones_like(x, dtype=int)

            # find the pixel coordinates for each unique value of the
            # region image
            ri = indices.reverse(r, ignore=(0,))
            for regid, pixid in ri.items():
                regid = int(regid)

                # pixels and weights for this region
                xx = x[pixid]
                yy = y[pixid]
                ww = w[pixid]

                # make and save the region
                reg = DispersedRegion(xx, yy, ww, self.segid, regid,
                                      ltv=self.ltv)
                self.append(reg)

            # parse the region image
            # ri = indices.reverse((seg == self.segid)*reg, ignore=(0,))
            # for regid, (yy, xx) in ri.items():
            #    ww = img[yy, xx]/self.norm
            #    reg = DispersedRegion(xx, yy, ww, self.segid, regid,
            #                          ltv=self.ltv)
            #     self.append(reg)

            #
            # with weights process if for different segmentation regions
            # segid < 0 ---> compound
            # segid > 0 ---> simple
            #
            # Eventually have "complex" morphologies in the spectral regions
            # that can describe an arbitrary arrangement of spectral regions
            #
            # if self.segid < 0:
            #     # a compound source, so each pixel in a separate DispersedRegion
            #     for i,(xx,yy,ww) in enumerate(zip(x,y,w),start=1):
            #         reg=DispersedRegion([xx],[yy],[ww],self.segid,i,
            #                           ltv=self.ltv)
            #         self.append(reg)
            # elif self.segid >0:
            #     self.append(DispersedRegion(x,y,w,self.segid,0,ltv=self.ltv))
            # else:
            #     LOGGER.warning(f"Ignoring {segid=}")

    def compute_weights(self, subimg, x, y, whttype='pixels', profile='gaussian',
                        negfunc='positivity', epsilon=1e-9):

        # set some variables
        self.profile = profile.lower()
        self.whttype = whttype.lower()

        # process the pixels based on the whttype specified
        if self.whttype == 'pixels':
            w = subimg[y, x]
        elif self.whttype == 'constant':
            w = np.ones_like(x, dtype=float)
        elif self.whttype == 'fitprofile':
            # get an estimate on the source's elliptical properties
            xc, yc, a, b, theta = self.ellipse_parameters(subimg, x, y)
            amplitude = np.amax(subimg[y, x])

            # define the fitting algorithm
            fitter = fitting.LevMarLSQFitter()

            # determine the spatial profile to fit
            if self.profile in ('gaussian', 'gauss'):
                mod = models.Gaussian2D(amplitude=amplitude,
                                        x_mean=xc, x_stddev=a,
                                        y_mean=yc, y_stddev=b,
                                        theta=theta)
            elif self.profile == 'sersic':
                mod = models.Sersic2D(amplitude=amplitude,
                                      x_0=xc, y_0=yc,
                                      n=2.5, r_eff=a,
                                      ellip=np.sqrt(1 - (b / a)**2),
                                      theta=theta)
            else:
                msg = f'profile {profile} is not supported'
                raise NotImplementedError(msg)

            # fit the data and compute the model
            p = fitter(mod, x, y, subimg[y, x])
            w = p(x, y)
        else:
            self.whttype = 'pixels'
            LOGGER.warning(f'WHTTYPE ({whttype}) not found, using {self.whttype}')
            w = subimg[y, x]

        # must remap the pixels into some positive weights
        self.negfunc = negfunc.lower()
        if self.negfunc == 'positivity':
            v = 0.5 * (w + np.sqrt(w * w + epsilon))
        elif self.negfunc == 'abs':
            v = np.abs(w)
        elif self.negfunc == 'relu':
            v = np.maximum(w, 0.)
        elif self.negfunc == 'minimum':
            v += np.amin(w)
        else:
            self.negfunc = 'positivity'
            LOGGER.warning(
                f'Negative pixel function ({negfunc.lower()}) not found, using {self.negfunc}')
            v = 0.5 * (w + np.sqrt(w * w + epsilon))

        # normalize
        v /= np.sum(v)

        return v

    def __str__(self):
        return super().__str__()

    @property
    def lamb0(self):
        """
        Property for the starting wavelength
        """
        return self.specpars['lamb0']

    @lamb0.setter
    def lamb0(self, l):
        """
        Setter for the starting wavelength

        Properties
        ----------
        l : int,float
            The starting wavelength
        """
        self.specpars['lamb0'] = l

    @property
    def lamb1(self):
        """
        Property for the ending wavelength
        """
        return self.specpars['lamb1']

    @lamb1.setter
    def lamb1(self, l):
        """
        Setter for the ending wavelength

        Properties
        ----------
        l : int,float
            The ending wavelength
        """
        self.specpars['lamb1'] = l

    @property
    def dlamb(self):
        """
        Property for the delta wavelength
        """
        return self.specpars['dlamb']

    @dlamb.setter
    def dlamb(self, l):
        """
        Setter for the delta wavelength

        Properties
        ----------
        l : int,float
            The delta wavelength
        """
        self.specpars['dlamb'] = l

    @property
    def is_compound(self):
        """
        A flag if the source is a compound source (ie. has multiple
        `DispersedRegion`s)
        """
        return len(self) > 1

    @property
    def pixelarea(self):
        """
        Compute the pixel area associated with this source
        """
        return wcsutils.proj_plane_pixel_area(self.wcs)

    @property
    def name(self):
        return str(self.segid)

    @property
    def nregions(self):
        """
        Compute the number of dispersed regions
        """
        return len(self)

    def pixels(self, **kwargs):  # weights=False,applyltv=False,dtype=int):
        """
        A generator to loop over all direct image pixels

        Parameters
        ----------
        applyltv : bool, optional
            Flag to apply the LTV before returning.  Internally to this
            source, the pixels are stored in lowest integers, and not on
            the original pixel grid.  Default is False

        dtype : type, optional
            The dtype to return the coordinates.  Default is int.

        Returns
        -------
        x : arb type
            The x-coordinates

        y : arb type
            The y-coordinates
        """

        for region in self:
            yield from region.pixels(**kwargs)

    def items(self):
        """
        Generator to implement `items()`

        Returns
        -------
        segid : int
             The segmentation ID

        regid : int
             The spectral region ID

        region : `DispersedRegion`
             The disersed region region
        """
        for i, region in enumerate(self):
            yield (self.segid, i), region

    def as_HDU(self, **kwargs):
        """
        Method to return a header-data unit (HDU) from this source

        Parameters
        ----------
        kwargs : dict, optional
            A set of keyword/value pairs to add to the header

        Returns
        -------
        hdu : `astropy.io.fits.ImageHDU`
            the output HDU

        Notes
        -----
        The EXTNAME keyword will be a string-version of the segid.
        If the source is compound, then it will be written with 3-dimensional
            WCS.
        """
        hdr = fits.Header()
        hdr['EXTNAME'] = (str(self.segid),)

        hdr['SEGID'] = (self.segid, 'Segmentation ID')
        hdr['GRPID'] = (kwargs.get('group', self.grpid), 'Group ID')
        hdr['RA'] = (self.adc[0], 'Right Ascension (deg)')
        hdr['DEC'] = (self.adc[1], 'Declination (deg)')
        hdr['X'] = (self.xyc[0], 'X barycenter')
        hdr['Y'] = (self.xyc[1], 'Y barycenter')
        hdr['FLUX'] = (self.flux, 'instrumental flux in direct image')
        hdr['MAG'] = (self.mag, 'AB magnitude in direct image')
        hdr['FNU'] = (self.fnu, 'flux in erg/s/cm2/Hz in direct image')
        hdr['NPIXELS'] = (self.npixels, 'total number of extracted pixels')
        hdr['WHTTYPE'] = (self.whttype, 'Type of source profile weights')
        if self.whttype == 'fitprofile':
            hdr['WHTPROF'] = (self.profile, 'The analytic profile for the weights')
        hdr['NREGIONS'] = (self.nregions, 'number of spectral regions')
        hdr['AREA'] = (self.npixels * self.pixelarea, 'source area (arcsec2)')
        headers.add_stanza(hdr, 'Source Properties', before='SEGID')

        hdr['BCKSUB'] = (self.local_back, 'Was local background in direct image subtracted')
        hdr['BCKSIZE'] = (self.backsize, 'Size of background annulus (in pix)')
        hdr['BCKVAL'] = (self.background, 'Background level in direct image')
        hdr['BCKLOSIG'] = (self.nsig[0], 'N sigma for low level')
        hdr['BCKHISIG'] = (self.nsig[1], 'N sigma for high level')
        headers.add_stanza(hdr, 'Direct Image Background', before='BCKSUB')

        if self.is_compound:
            hdr['NAXIS'] = (3, 'number of WCS axes')
            hdr['NAXIS1'] = self.wcs._naxis[0]
            hdr['NAXIS2'] = self.wcs._naxis[1]
            # hdr['NAXIS3']=nlam
            raise RuntimeError("gotta get nlam")

            hdr['CTYPE1'] = self.wcs.ctype[0]
            hdr['CTYPE2'] = self.wcs.ctype[1]
            hdr['CTYPE3'] = 'wavelength (A)'

            hdr['CRPIX1'] = self.wcs.crpix[0]
            hdr['CRPIX2'] = self.wcs.crpix[1]
            hdr['CRPIX3'] = 1

            hdr['CRVAL1'] = self.wcs.crval[0]
            hdr['CRVAL2'] = self.wcs.crval[1]
            # hdr['CRVAL3']=lam0

            hdr['LTV1'] = self.ltv[0]
            hdr['LTV2'] = self.ltv[1]

            hdr['CD1_1'] = self.wcs.cd[0, 0]
            hdr['CD1_2'] = self.wcs.cd[0, 1]
            hdr['CD2_1'] = self.wcs.cd[1, 0]
            hdr['CD2_2'] = self.wcs.cd[1, 1]
            # hdr['CD3_3']=dlam

            hdr['CUNIT1'] = self.wcs.cunit[0]
            hdr['CUNIT2'] = self.wcs.cunit[1]
            hdr['CUNIT3'] = 'A'

            hdr['DISPAXIS'] = (3, 'dispersion axis')
            headers.add_stanza(hdr, 'WORLD-COORDINATE SYSTEM', before='NAXIS')

            # set the data
            dat = np.zeros((hdr['NAXIS3'], hdr['NAXIS2'], hdr['NAXIS1']),
                           dtype=float)

            hdu = fits.ImageHDU(dat, hdr)
        else:
            # set the data as a single spectrum
            dat = self[0].sed.data
            hdu = fits.BinTableHDU(dat, hdr)
        return hdu

    def image_coordinates(self, x, y, dtype=None):
        """
        Method to transform coordinates by the LTV keywords

        Parameters
        ----------
        x : int, float, or `np.ndarray`
            The x-coordinates

        y : int, float, or `np.ndarray`
            The y-coordinates

        dtype : type or None, optional
            Variable to recast the data types.  If None, then no retyping
            is done.  Default is None

        Returns
        -------
        xx : arbitrary type
            The transformed x-coordinates

        yy : arbitrary type
            The transformed y-coordinates

        Notes
        -----
        The input coordinates must have the same shape, their dtype does
        not matter.  The output coordinates will have that shape.
        """

        xx = x - self.ltv[0]
        yy = y - self.ltv[1]
        if dtype is not None:
            xx = xx.astype(dtype)
            yy = yy.astype(dtype)

        return xx, yy

    def set_spectral_parameters(self, **kwargs):
        """
        Method to set the spectral parameters for all of the spectral regions
        that are part of this source

        Parameters
        ----------
        kwargs : dict, optional
            Dictionary of keywords, can be 'wave0','wave1','dwave','scale'

        """

        for region in self:
            region.set_spectral_parameters(**kwargs)

    def get_spectral_parameters(self):
        """
        Method to get the spectral parameters from this object's attributes

        Returns
        -------
        pars : list
            A list of tuples of the extraction parameters
        """
        raise NotImplementedError()

        # pars=[region.get_spectral_parameters() for region in self]
        #
        # dtypes=[for p in pars[0]]
        #
        # return pars

    def load_sedlib(self, sedlib, throughput=None):
        """
        Method to load SEDs from an SED library

        Parameters
        ----------
        sedlib : `SEDFile`
            The library of SEDs

        throughput : `su.core.photometry.Throughput` or None
            The throughput curve used to normalized the SED.  If set to
            `None`, then do not normalize.

        """
        for regkey, region in self.items():
            sed = sedlib[regkey]

            if sed and throughput:
                sed.normalize(throughput, self.fnu * np.sum(region.w))

            region.sed = sed

            sed.write_file(f'spectra/{self.segid}_{region.regid}.csv')

    def load_sed_images(self):
        pass

    def load_sed_array(self, waves, flam):
        """
        Method to set all the SEDs of all the spectral regions to this
        spectrum, but scaled by their relative flux levels

        Parameters
        ----------
        waves : `np.ndarray`
            Wavelength array

        flam : `np.ndarray`
            flux array in flam units

        Notes
        -----
        The input spectrum is used for all the `DispersedRegion`s
        """

        for regkey, region in self.items():
            region.sed.append(waves, flam * np.sum(region.w))

    def write_seds(self, filetype='sed', path=None, **kwargs):
        """
        Method to write the SEDs to disk

        Parameters
        ----------
        filetype : str, optional
            The type of file to write.  See `su.core.photometry.SED` for the
            possible file types.  Default is 'sed'

        path : str, optional
            The path to write the SEDs to.  Will create the output directory
            if it is set, but does not exist.  Default is '.'

        kwargs : dict, optional
            Additional keyword/value pairs that are passed into
            `su.core.photometry.SED().write_file`

        """

        if path is None:
            path = '.'
        else:
            if not os.path.isdir(path):
                os.mkdir(path)

        for regid, region in enumerate(self):
            filename = os.path.join(path,
                                    f'segid_{self.segid}_{regid}.{filetype}')
            region.sed.write_file(filename, **kwargs)

    @staticmethod
    def ellipse_parameters(img, x, y):
        """
        Staticmethod to compute the ellipse parameters from an image

        Parameters
        ----------
        img : `np.ndarray`
            A 2d array that describes the pixel brightnesses of the source

        wht : `np.ndarray`
            A 2d array that describes the weights.  In the simplest case, this
            can be a boolean array, such as from a segmentation map, but could
            also include uncertainties.

        Returns
        -------
        x_mean : float
            Expected X-centroid (in pixels)

        y_mean : float
            Expected Y-centroid (in pixels)

        a : float
            semimajor axis (in pixels)

        b : float
            semiminor axis (in pixels)

        theta : float
            orientation angle (in radians)

        Notes
        -----
        This follows the definitions in the Source Extractor manual

        """

        weights = img[y, x]

        # compute various moments
        X = np.average(x, weights=weights)
        Y = np.average(y, weights=weights)
        X2 = np.average(x * x, weights=weights) - X * X
        Y2 = np.average(y * y, weights=weights) - Y * Y
        XY = np.average(x * y, weights=weights) - X * Y

        # compute the difference and average of the 2nd order
        dif = (X2 - Y2) / 2.
        ave = (X2 + Y2) / 2.

        # compute the main radical
        rad = np.sqrt(dif**2 + XY**2)

        # compute axes
        a = np.sqrt(ave + rad)
        b = np.sqrt(ave - rad)

        # compute the position angle (in radians)
        theta = 0.5 * np.arctan2(XY, dif)

        return X, Y, a, b, theta

    def plot(self, minfactor=1e-3, padding=1):
        """
        Method to make a plot.

        Parameters
        ----------
        minfactor : float, optional
            The minimum display value (in units of peak weight).   Default is 1e-3

        padding : int, optional
            Amount to pad the images.  Default is 1.
        """

        fig, axes = plt.subplots(1, 2, sharex=True, sharey=True)

        x = []
        y = []
        regid = []
        weight = []
        for i, region in enumerate(self, start=1):
            x.extend(region.x)
            y.extend(region.y)
            regid.extend([i] * len(x))
            weight.extend(region.w)
        x = np.asarray(x, dtype=int)
        y = np.asarray(y, dtype=int)
        x0, x1 = np.amin(x) - padding, np.amax(x) + padding
        y0, y1 = np.amin(y) - padding, np.amax(y) + padding
        nx, ny = x1 - x0 + 1, y1 - y0 + 1
        reg = np.full((ny, nx), np.nan)
        img = np.full((ny, nx), np.nan)

        # reg[y - y0 + padding, x - x0 + padding] = regid
        # img[y - y0 + padding, x - x0 + padding] = weight
        reg[y - y0, x - x0] = regid
        img[y - y0, x - x0] = weight

        mx = np.nanmax(img)
        axes[0].imshow(reg, origin='lower')
        axes[1].imshow(img, origin='lower', norm=colors.LogNorm(mx * minfactor, mx))

        plt.tight_layout()
        plt.show()
