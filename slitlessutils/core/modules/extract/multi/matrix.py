from collections import namedtuple

import numpy as np
from scipy import sparse
from tqdm import tqdm

from .....config import Config
from .....logger import LOGGER
from ....tables import PDTFile
from ....utilities import as_iterable, headers  # , indices
from .lcurve import LCurve
from .result import Result


ImageData = namedtuple("ImageData", ['dataset', 'detname', 'pixels'])


class Matrix:

    def __init__(self, extorders, **kwargs):
        # extract the defaults
        mskorders = kwargs.get('mskorders', None)
        invmethod = kwargs.get('invmethod', 'lsqr')
        path = kwargs.get('path', 'su_tables')
        minunc = kwargs.get('minunc', 1e-10)

        LOGGER.debug("must finish documentation")

        # save some things
        self.path = path
        self.invmethod = invmethod
        self.minunc = minunc
        self.extorders = as_iterable(extorders)
        self.mskorders = as_iterable(mskorders)
        self.fluxscale = Config().fluxscale

    def __len__(self):
        return self.A.A.nnz if self else 0

    def __bool__(self):
        return hasattr(self, 'A') and (self.A.A.nnz > 0)

    def __str__(self):
        return f'Sparse matrix with {len(self)} elements'

    def load_image_matrix(self, h5, detdata, sources, gpx):
        '''
        Load a CSR matrix for a given WFSS image

        Inputs
        ------
        h5 : h5py file object
            The file to load the table from

        detdata : `WFSSDetector`
            The detector to load the table for

        sources : `SourceCollection`
            The source collection for this image

        gpx : `np.ndarray`
            A mask of good pixels

        Returns
        -------
        A : `scipy.sparse.csr_matrix`
            The CSR sparse matrix for this image

        yx : tuple of `np.ndarray`s
            The pixel coordinates in the WFSS image for which we have data

        '''

        # will need this
        npix = detdata.npixels()

        # get the flat field
        flatfield = detdata.config.load_flatfield()

        # a container to hold all the submatrices
        A = []
        for ordname in self.extorders:
            h5.load_order(ordname)
            sensitivity = detdata.config[ordname].sensitivity

            for source in sources.values():

                # get the extration properties for this source
                if hasattr(source, 'extpars'):
                    extpars = source.extpars
                else:
                    extpars = self.defpars

                # fetch the limits and the number of elements
                limits = extpars.limits()
                nlam = len(extpars)

                # the dimensionality of this local submatrix
                shape = (npix, nlam)

                for regidx, region in enumerate(source):
                    tab = h5.load_odt(region)
                    if tab:
                        wav = tab.wavelengths()
                        x = tab.get('x', dtype=int)
                        y = tab.get('y', dtype=int)
                        val = tab.get('val')

                        # remove badpixels
                        g = np.where(gpx[y, x])[0]
                        if g.size > 0:
                            wav = wav[g]
                            x = x[g]
                            y = y[g]
                            val = val[g]

                            # matrix coordinates
                            i = np.ravel_multi_index((y, x), dims=detdata.shape)
                            j = np.digitize(wav, limits) - 1

                            # compute weights
                            sens = sensitivity(wav) * self.fluxscale
                            flat = flatfield(x, y, wav)
                            area = detdata.relative_pixelarea(x, y)

                            # apply calibs to the weights
                            aij = (val * sens * flat * area)

                            # make the local sparse matrix
                            Asub = sparse.csr_matrix((aij, (i, j)), shape=shape)

                        else:
                            # create an empty matrix in case all pixels got
                            # masked out for whatever reason
                            Asub = sparse.csr_matrix(shape)

                    else:
                        # create an empty matrix for objects out of the field
                        Asub = sparse.csr_matrix(shape)

                    # save it to the list, so we can stack them later
                    A.append(Asub)

        # stack all the sparse matrices
        A = sparse.hstack(A)

        # find the rows that are all zero
        yx = np.where(np.diff(A.indptr))[0]

        # remove the zero rows
        A = A[yx]

        # get the image coordinates
        y, x = np.unravel_index(yx, shape=detdata.shape)

        return A, (y, x)

    def build_matrix(self, data, sources, group=0):
        '''
        Build a sparse matrix and ancillary data

        Inputs
        ------
        data : `WFSSCollection`
            The WFSS data

        sources : `SourceCollection`
            The sources to extract

        grpid : int
            The group that these sources refer to.  Default=0

        Returns
        -------
        None

        '''

        # do some quick checks
        if not self.extorders:
            LOGGER.warning('no valid extraction orders')
            return

        # we code it this way, so we can instantiate a matrix object
        # that contains all the defaults to update_headers.
        self.group = group

        # get some dimensionalities
        self.nfiles = len(data)
        self.nsources = len(sources)
        self.segids = list(sources.keys())
        self.defpars = data.get_parameters()

        # issue a message
        LOGGER.info("Building matrix with\n" +
                    f"       WFSS files: {self.nfiles}\n" +
                    f"          Sources: {self.nsources}\n" +
                    f"         Group ID: {self.group}")

        # lists for output data
        # b = A.x, so b, and A are the linear variables
        # pixels is a list of all the image pixels we've swept up
        self.A = []
        self.b = []
        self.imgdata = []

        # process each WFSS image
        for datum in tqdm(data, desc="Loading WFSS Images", total=self.nfiles,
                          dynamic_ncols=True):

            # load the data for this WFSS image
            with PDTFile(datum, path=self.path, mode='r') as h5:

                for detindex, (detname, detdata) in enumerate(datum.items()):
                    h5.load_detector(detname)

                    # read the images
                    sci, hdr = detdata.readfits('science', header=True)
                    unc = detdata.readfits('uncertainty')
                    dqa = detdata.readfits('dataquality')

                    # adjust the units if necessary
                    bunit = hdr.get('BUNIT', 'electron/s')
                    if bunit.lower() in ('electron', 'electrons', 'e', 'e-'):
                        phdr = detdata.primaryheader()
                        time = phdr['EXPTIME']
                        sci /= time
                        unc /= time

                    # initialize a good-pixel mask as the pixels that are finite
                    gpx = np.isfinite(sci) & np.isfinite(unc)

                    # update the good-pixel mask with a bitmask
                    if detdata.config.bitmask:
                        gpx &= (np.bitwise_and(dqa, detdata.config.bitmask) == 0)

                    # if detdata.config.bitmask:
                    #    gpx = np.bitwise_and(dqa, detdata.config.bitmask) == 0
                    # else:
                    #    gpx = np.ones_like(dqa, dtype=bool)

                    # this will work, but maybe less efficient?
                    # bad = np.logical_not(np.isfinite(sci) & np.isfinite(unc))
                    # gpx[bad] = False

                    # update the mask for the mask orders
                    for ordname in self.mskorders:
                        h5.load_order(ordname)
                        for source in sources.values():
                            omt = h5.load_omt(source)
                            gpx[omt.get('y'), omt.get('x')] = False

                    # load all the submatrices for each order
                    Aimg, yx = self.load_image_matrix(h5, detdata, sources, gpx)

                    # grab the pixel values from data and uncertainty
                    # but NB: the unc could potentially be zero, which will
                    # cause problems with weighting by inverse-variance.
                    # so the value minunc will set the minimum uncertainty
                    # to avoid this problem
                    s = sci[yx]
                    u = np.maximum(unc[yx], self.minunc)

                    # normalize the data and matrices by the uncertainties
                    # for inverse-variance weighting of the solution (ie.
                    # solving maximum likelihood).  But will change to a
                    # column matrix for latter use
                    Aimg /= u[:, np.newaxis]
                    s /= u

                    # save the data and matrices
                    self.b.extend(s)
                    self.A.append(Aimg.tocsc())

                    # save some light-weight things about this image
                    imgdata = ImageData(datum.dataset, detname, yx)
                    self.imgdata.append(imgdata)

        # store things in a way to use later
        self.A = sparse.linalg.aslinearoperator(sparse.vstack(self.A))
        self.b = np.asarray(self.b)
        self.frob = sparse.linalg.norm(self.A.A)

        # make a mask for spectral elements that have no constraints in the
        # data, which corresponds to columns that are all zero
        self.unconstrained = (np.diff(self.A.A.indptr) == 0)

        # set some variables
        self.npix, self.npar = self.A.shape
        nel = self.npix * self.npar
        if nel == 0:
            self.density = np.nan
        else:
            self.density = float(self.A.A.nnz) / nel

        # build an L-Curve object to store the data
        self.lcurve = LCurve(norm=self.frob)

        if self:
            # build a target spectrum
            target = []
            for segid, source in sources.items():
                if hasattr(source, 'extpars'):
                    waves = source.extpars.wavelengths()
                else:
                    waves = self.defpars.wavelengths()

                for regid, region in source.items():
                    target.extend(region.sed(waves))

            target = np.asarray(target)
            if len(target) != self.npar:
                LOGGER.error("Invalid target spectrum")
                exit()

            # set the damping target
            self.set_damping_target(target)
        else:
            LOGGER.warning('The matrix is empty, no calculations possible.')

        LOGGER.info("Finished building matrix")

    def compute_model(self, x):
        """
        Method to compute the model, which is given as a matrix-product
        between this object and a vector of spectral values

        Parameters
        ----------
        x : `np.ndarray`
           An array of spectra that have been concatenated in the same way
           that this `Matrix` was created.  This is used to estimate the
           expected signal on the detector(s) given a known spectrum

        Returns
        -------
        out : `np.ndarray`
           The expected signal given by

           .. math::
              b = A\\dot x

           In matrix notation this is
           .. math::
              b_i = \\sum_j A_{i,j}x_{j}

           Hence the name of the variable `xj`.

        Notes
        -----
        This allows for quick simulations or a chi2-like test, given a
        set of spectra

        """

        mod = {}
        if self:
            b = self.A.matvec(x)

            s = slice(0, 0)
            for imgdata in self.imgdata:
                s = slice(s.stop, s.stop + len(imgdata.pixels[0]))

                key = (imgdata.dataset, imgdata.detname)
                mod[key] = list(zip(imgdata.pixels[0], imgdata.pixels[1], b[s]))

        return mod

    def set_damping_target(self, target):
        """
        Method to set the damping target

        Parameters
        ----------
        target : `np.ndarray`
           The damping target as a concatenated array of spectra

        """
        if self:
            self.target = target / self.fluxscale
            self.b -= self.A.matvec(self.target)
        else:
            LOGGER.warning("Cannot set damping target until matrix is built")

    # --------------------------------------------------------------------
    #
    #              functions for solving the linear problem
    #
    # --------------------------------------------------------------------

    def invert(self, logdamp, scale=True, **kwargs):
        """
        Method to the matrix inversion

        Parameters
        ----------
        logdamp : float or None
           The log10 of the damping parameter.  If None, then damp=0.0
           is understood

        scale : bool, optional
           A flag that the damping parameter should be rescaled by the
           Frobenius norm of the matrix.  See Ryan, Casertano, Pirzkal 2018
           https://ui.adsabs.harvard.edu/abs/2018PASP..130c4501R/abstract
           for discussion of this choice.

        kwargs : dict, optional
           Additional keywords passed to either
           `scipy.sparse.linalg.lsqr()` or `scipy.sparse.linalg.lsmr()`

        """

        if logdamp is None:
            damp = 0.
        else:
            damp = 10.**logdamp
            if scale:  # and self.frob >1:
                damp *= self.frob

        LOGGER.info(f"Running {self.invmethod} with \u2113={damp}")

        # here calls the function
        r = self._invfunc(damp, **kwargs)

        # take the frobenius norm and targets out?
        if scale:  # and self.frob>1:
            r.damp /= self.frob

        if hasattr(self, 'target'):
            r += self.target

        # update the Lcurve
        self.lcurve.append(r.r1norm, r.xnorm, r.logdamp)

        return r

    def _lsqr(self, damp, **kwargs):
        """
        Helper function to call LSQR

        Parameters
        ----------
        damp : float
            The scaled damping target to use

        kwargs : dict, optional
            Additional arguments passed to `scipy.sparse.linalg.lsqr()`

        Returns
        -------
        r : `Result`
            The packaged result

        Notes
        -----
        Although possible, it's not generally necessary to directly call
        this function, but rather go through `self.invert()`

        """

        r = sparse.linalg.lsqr(self.A, self.b, damp=damp, calc_var=True, **kwargs)
        (x, istop, itn, r1norm, r2norm, anorm, acond, arnorm, xnorm, var) = r

        r = Result('LSQR', x, istop, itn, r1norm, r2norm, anorm, acond,
                   arnorm, xnorm, np.sqrt(var), damp)

        return r

    def _lsmr(self, damp, **kwargs):
        """
        Helper function to call LSMR

        Parameters
        ----------
        damp : float
            The scaled damping target to use

        kwargs : dict, optional
            Additional arguments passed to `scipy.sparse.linalg.lsmr()`

        Returns
        -------
        r : `Result`
            The packaged result

        Notes
        -----
        Although possible, it's not generally necessary to directly call
        this function, but rather go through `self.invert()`

        """

        r = sparse.linalg.lsmr(self.A, self.b, damp=damp, **kwargs)
        (x, istop, itn, norm, arnorm, anorm, acond, xnorm) = r

        r1norm = arnorm
        r2norm = np.sqrt(arnorm**2 + (damp * xnorm)**2)
        unc = np.full_like(x, np.nan)
        r = Result('LSMR', x, istop, itn, r1norm, r2norm, anorm, acond,
                   arnorm, xnorm, unc, damp)

        return r

    def compute_uncertainty(self):
        ones = np.ones_like(self.b, dtype=float)
        A2 = self.A.A.copy()
        A2.data = np.square(A2.data)

        r = sparse.linalg.lsqr(A2, ones)
        return np.sqrt(r[0])

    # --------------------------------------------------------------------
    #
    #                      some properites
    #
    # --------------------------------------------------------------------

    @property
    def is_underdetermined(self):
        return self.npar > self.npix

    @property
    def is_determined(self):
        return self.npar == self.npix

    @property
    def is_overdetermined(self):
        return self.npar < self.npix

    @property
    def invmethod(self):
        return self._invmethod

    @invmethod.setter
    def invmethod(self, invmethod):
        im = invmethod.lower()
        if im == 'lsqr':
            self._invfunc = self._lsqr
            self._invmethod = im
        elif im == 'lsmr':
            self._invfunc = self._lsmr
            self._invmethod = im
        else:
            LOGGER.warning(f'{invmethod=} is not supported. Using `LSQR`')
            self._invfunc = self._lsqr
            self._invmethod = 'lsqr'

    # --------------------------------------------------------------------
    #
    #                      outputting functions
    #
    # --------------------------------------------------------------------

    def update_header(self, hdr):
        """
        Method to update a fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            The header to update
        """

        if not self:
            LOGGER.warning("Cannot update header, matrix uninitialized")
            return

        hdr['NNZ'] = (len(self), 'number of non-zero matrix elements')
        hdr['NPIX'] = (self.npix, 'number of pixels analyzed (knowns)')
        hdr['NPAR'] = (self.npar, 'number of parameters measured (unknowns)')
        hdr['NDOF'] = (self.npix - self.npar, 'number of degrees of freedom')
        hdr['DENSITY'] = (self.density, 'fraction of non-zero elements')
        hdr['FROBNORM'] = (self.frob, 'Frobenius norm')

        hdr['INVMETH'] = (self.invmethod, 'Method to compute inverse')
        hdr['MINUNC'] = (self.minunc, 'minimum uncertainty in e/s')
        hdr['EXTORDER'] = (','.join(self.extorders), 'extraction orders')
        hdr['MSKORDER'] = (','.join(self.mskorders), 'masked orders')
        headers.add_stanza(hdr, 'Matrix Properties', before='NNZ')
