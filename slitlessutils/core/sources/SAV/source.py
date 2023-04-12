import os

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from astropy.wcs import WCS, utils as wcsutils
import matplotlib.pyplot as plt
import numpy as np
from skimage.segmentation import expand_labels

from ...logger import LOGGER
from .spectralregion import SpectralRegion
from ..utilities import headers, indices


class Source(list):
    """
    Class to hold a single spectral source, which in turn may contain any
    number of spectral regions that have a unique spectrum.

    inherits from list.  Each element of the list is a `SpectralRegion`,
    which controls the unique spectrum.
    """

    # maximum number of spectral regions
    MAXSPECREG = 1000

    def __init__(self, img, seg, hdr, local_back=True, backsize=5, nsig=(5, 5),
                 whttype='pixels', profile='gaussian', zeropoint=26., grpid=0,
                 specregions=None):
        """
        Initializer

        Parameters
        ----------
        img : `np.ndarray`
            A two-dimensional array from the direct image

        seg : `np.ndarray`
            A two-dimensional array from the segmentation image

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
            is a tuple, then the numbers are interepreted as (low,high)
            values. If it is a single value, then will be both low and high
            clipping values.  Default is (5.,5.).

        whttype : str, optional
            Way to control how the direct image values are processed into
            the weights.  If equal to:

            'pixels'       Use the absolute value of the pixels

        profile : str, optional
            Name of the profile to fit.  See `astropy.modeling.models`.
            This is only used if ```whttype='fitprofile'```.
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

        specregions : `np.ndarray` or None, optional
            If an np.array, then this is a "segmentation map" like object,
            then these pixels get aggregated.


        Notes
        ----
        This is a cornerstone object to slitlessutils.

        """

        self.wcs = WCS(hdr, relax=True)
        self.ltv = [hdr.get('LTV1', 0.), hdr.get('LTV2', 0.)]

        # record some things
        self.segid = hdr['SEGID']
        self.grpid = grpid
        self.whttype = whttype
        self.backsize = backsize
        self.local_back = local_back

        # set the spectral parameters
        self.specpars = {'lamb0': None, 'lamb1': None, 'dlamb': None}

        # test the number of sigma for the sky
        if isinstance(nsig, (float, int)):
            nsig = (nsig, nsig)
        self.nsig = nsig

        # find the pixels for this segmentation
        y, x = np.where((img > 0) & (seg == self.segid))
        self.npixels = len(x)
        if self.npixels > 0:
            # compute and subtract the local sky background:
            if self.local_back and self.backsize > 0:
                # get a bounding box
                x0 = max(np.amin(x)-self.backsize, 0)
                y0 = max(np.amin(y)-self.backsize, 0)
                x1 = min(np.amax(x)+self.backsize, img.shape[1]-1)
                y1 = min(np.amax(y)+self.backsize, img.shape[0]-1)

                # cut out regions
                subimg = img[y0:y1, x0:x1]
                subseg = seg[y0:y1, x0:x1]

                # grow the region to find the local sky background
                # but use a mask in case it results in empty array
                subseg2 = expand_labels(subseg, distance=self.backsize)
                mask = np.logical_or(subseg2 != self.segid, subseg != 0)
                try:
                    ave, med, sig = sigma_clipped_stats(subimg, mask=mask,
                                                        sigma_lower=self.nsig[0],
                                                        sigma_upper=self.nsig[1])
                except BaseException:
                    ave = 0.0

                self.background = ave
            else:
                self.background = 0.0

            img -= self.background

            # process the pixels based on the whttype specified
            if self.whttype == 'pixels':
                w = np.abs(img[y, x])

            elif self.whttype == 'fitprofile':

                # get an estimate on the source's elliptical properties
                subwht = (subseg == self.segid).astype(float)
                xc, yc, a, b, theta = self.ellipse_parameters(subimg, subwht)
                amplitude = np.amax(subimg[y, x])

                # define the fitting algorithm
                fitter = fitting.LevMarLSQFitter()

                # determine the spatial profile to fit
                self.profile = profile.lower()
                if self.profile == 'gaussian':
                    p_init = models.Gaussian2D(amplitude=amplitude,
                                               x_mean=xc, x_stddev=a,
                                               y_mean=yc, y_stddev=b,
                                               theta=theta)
                elif self.profile == 'sersic':
                    p_init = models.Sersic2D(amplitude=amplitude,
                                             x_0=xc, y_0=yc,
                                             n=2.5, r_eff=a,
                                             ellip=np.sqrt(1-(b/a)**2),
                                             theta=theta)

                else:
                    msg = f'profile {profile} is not supported'
                    raise NotImplementedError(msg)

                # fit the data and compute the model
                p = fitter(p_init, x, y, img[y, x])
                w = p(x, y)

            else:
                self.whttype = 'pixels'
                w = np.abs(img[y, x])

            # normalize the weights
            self.norm = np.sum(w)
            w /= self.norm

            # compute the centeroids
            xyc = (np.average(x, weights=w), np.average(y, weights=w))
            self.xyc = self.image_coordinates(*xyc)
            adc = self.wcs.all_pix2world(*xyc, 0)
            self.adc = (float(adc[0]), float(adc[1]))

            # compute the area of the object
            self.area = self.npixels*self.pixelarea

            # compute some fluxes
            self.flux = np.sum(img[y, x])

            self.fnu = self.flux*10.**(-0.4*(zeropoint+48.6))
            if self.flux < 0:
                self.mag = 99.
            else:
                self.mag = -2.5*np.log10(self.flux)+zeropoint

            # put coordinates back on original footprint
            # x,y=self.image_coordinates(x,y,dtype=np.int)

            #
            # with weights process if for different segmentation regions
            # segid < 0 ---> compound
            # segid > 0 ---> simple
            #
            # Eventually have "complex" morphologies in the spectral regions
            # that can describe an arbitrary arrangement of spectral regions
            #
            if self.segid < 0:
                # a compound source, so each pixel in a separate SpectralRegion
                for i, (xx, yy, ww) in enumerate(zip(x, y, w), start=1):
                    reg = SpectralRegion([xx], [yy], [ww], self.segid, i, ltv=self.ltv)
                    self.append(reg)
            elif self.segid > 0:
                self.append(SpectralRegion(x, y, w, self.segid, 0, ltv=self.ltv))
            else:
                LOGGER.warning(f"Ignoring {segid=}")

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
        `SpectralRegions`)
        """
        # return self.nregions >1
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
        Compute the number of spectral regions
        """
        # return len(self.spectralregions)
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

        # if applyltv:
        #    for region in self:
        #        for x,y in region.pixels():
        #            yield self.image_coordinates(x,y,dtype=dtype)
        # else:
        #    for region in self:
        #        for x,y in region.pixels():
        #            yield (dtype(x),dtype(y))

        # if applyltv:
        #    #for region in self.spectralregions:
        #    for region in self:
        #        for x,y in region.pixels():
        #            #if applyltv:
        #            yield self.image_coordinates(x,y,dtype=dtype)
        # else:
        #    #for region in self.spectralregions:
        #    for region in self:
        #        yield from region.pixels()

    def items(self):
        """
        Generator to implement `items()`

        Returns
        -------
        segid : int
             The segmentation ID

        regid : int
             The spectral region ID

        region : `SpectralRegion`
             The spectral region
        """

        # for i,region in enumerate(self.spectralregions):
        #    yield (self.segid,i),region

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
        The EXTNAME kwyord will be a string-version of the segid.
        If the source is compound, then it wil be written with 3-dimensional
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
        hdr['AREA'] = (self.npixels*self.pixelarea, 'source area (arcsec2)')
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
            # dat=self.spectralregions[0].sed.data
            dat = self[0].sed.data
            hdu = fits.BinTableHDU(dat, hdr)
        return hdu


#    def update_header(self,hdr,group=0):
#
#        hdr['SEGID']=(self.segid,'Segmentation ID')
#        hdr['GRPID']=(group,'Group ID')
#        hdr['RA']=(self.adc[0],'Right Ascension (deg)')
#        hdr['DEC']=(self.adc[1],'Declination (deg)')
#        hdr['X']=(self.xyc[0],'X barycenter')
#        hdr['Y']=(self.xyc[1],'Y barycenter')
#        hdr['FLUX']=(self.flux,'instrumental flux in direct image')
#        hdr['MAG']=(self.mag,'AB magnitude in direct image')
#        hdr['FNU']=(self.fnu,'flux in erg/s/cm2/Hz in direct image')
#        hdr['NPIXELS']=(self.npixels,'total number of extracted pixels')
#        hdr['WHTTYPE']=(self.whttype,'Type of source profile weights')
#        hdr['NREGIONS']=(self.nregions,'number of spectral regions')
#        hdr['AREA']=(self.npixels*self.pixelarea,'source area (arcsec2)')
#        headers.add_stanza(hdr,'Source Properties',before='SEGID')
#
#
#        hdr['BCKSUB']=(self.local_back,'Was local background in direct image su#btracted')
#        hdr['BCKSIZE']=(self.backsize,'Size of background annulus (in pix)')
#        hdr['BCKVAL']=(self.background,'Background level in direct image')
#        hdr['BCKLOSIG']=(self.nsig[0],'N sigma for low level')
#        hdr['BCKHISIG']=(self.nsig[1],'N sigma for high level')
#        headers.add_stanza(hdr,'Direct Image Background',before='BCKSUB')
#
#
#
#        #for k,v in kwargs.items():
#        #    hdr[k]=v
#
#


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

        xx = x-self.ltv[0]
        yy = y-self.ltv[1]
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
                sed.normalize(throughput, self.fnu*np.sum(region.w))

            region.sed = sed

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
        The input spectrum is used for all the `SpectralRegion`s
        """

        for regkey, region in self.items():

            # region.sed.set_sed(waves,flam*np.sum(region.w))
            region.sed.append(waves, flam*np.sum(region.w))

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

        # for regid,region in enumerate(self.spectralregions):
        for regid, region in enumerate(self):
            filename = os.path.join(path, f'segid_{self.segid}_{regid}.{filetype}')
            region.sed.write_file(filename, **kwargs)

    @staticmethod
    def ellipse_parameters(img, wht):
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

        # find the good points and normalize the weights
        y, x = np.where(wht > 0)
        w = img[y, x]*wht[y, x]
        w /= np.sum(w)

        # compute the various moments
        X = np.sum(x*w)
        Y = np.sum(y*w)
        X2 = np.sum(x*x*w)-X*X
        Y2 = np.sum(y*y*w)-Y*Y
        XY = np.sum(x*y*w)-X*Y

        # compute the difference and average of the 2nd order
        dif = (X2-Y2)/2.
        ave = (X2+Y2)/2.

        # compute the main radical
        rad = np.sqrt(dif**2+XY**2)

        # compute axes
        a = np.sqrt(ave+rad)
        b = np.sqrt(ave-rad)

        # compute the position angle (in radians)
        theta = 0.5*np.arctan2(XY, dif)

        return X, Y, a, b, theta

    def plot(self):
        """
        Method to make a plot.

        This is not finished and should not be used
        """

        raise NotImplementedError

        x = []
        y = []
        w = []
        r = []
        # for regid,region in enumerate(self.spectralregions):
        for regid, region in enumerate(self):
            x.extend(region.x)
            y.extend(region.y)
            w.extend(region.w)
            r.extend([regid]*len(region))

        x = np.array(x)
        y = np.array(y)
        w = np.array(w)
        r = np.array(r)

        x0, x1 = np.amin(x), np.amax(x)
        y0, y1 = np.amin(y), np.amax(y)
        shape = (y1-y0+1, x1-x0+1)

        wht = np.full(shape, np.nan)
        reg = np.full(shape, np.nan)
        wht[y, x] = w
        reg[y, x] = r

        #
        #
        # x0,x1=np.inf,0
        # y0,y1=np.inf,0
        #
        # for region in self.spectralregions:
        #     x0=min(x0,np.amin(region.x))
        #     x1=max(x1,np.amax(region.x))
        #     y0=min(y0,np.amin(region.y))
        #     y1=max(y1,np.amax(region.y))
        #
        # shape=(y1-y0+1,x1-x0+1)
        # wht=np.full(shape,np.nan,dtype=float)
        # reg=np.full(shape,np.nan,dtype=float)
        #
        # for regid,region in enumerate(self.spectralregions):
        #     wht[region.y,region.x]=region.w
        #     reg[region.y,region.x]=regid

        fig, axes = plt.subplots(1, 2, sharex=True, sharey=True)
        axes[0].imshow(wht, origin='lower', interpolation='nearest')
        axes[0].set_title('weights')
        axes[1].imshow(reg, origin='lower', interpolation='nearest')
        axes[1].set_title('region')
        plt.tight_layout()
        fig.canvas.draw()
        plt.show()
