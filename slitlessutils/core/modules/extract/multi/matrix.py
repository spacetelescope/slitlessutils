import numpy as np
from scipy.sparse import coo_matrix, linalg
from tqdm import tqdm

from .....config import Config
from .....logger import LOGGER
from ....tables import PDTFile
from ....utilities import as_iterable, headers, indices
from .lcurve import LCurve
from .result import Result

# Just for notation sake
# variables with a 'g' are in the grism image
#                a 'u' are unique versions (so duplicates have been removed
#                      or summed over as applicable)
#                'comp' are compressed over the indices
#                'uniq' are uniq such that this is true: i=iuniq[icomp].
#                       which means that max(icomp)=len(iuniq)
# pixel coordinates will be represnted as either 1-d pixels (xyg) or
# a 2-d pixel pair (xg,yg).  This is because 1d coordinates are easier to
# work with in terms of summations, but harder when retrieving pixel
# values out of an image


class Matrix:
    """
    Class to implement the linear operator

    Parameters
    ----------
    extorders : str or list
       The spectral orders to extract.  If a list, tuple or set, then
       multiple orders will aggregated.

    mskorders : str or list or None, optional
       The spectral orders to mask in the extraction.  If a list, tuple,
       or set, then the multiple orders will be masked.

    invmethod : str, optional
       The inversion algorithm, options are 'lsqr' or 'lsmr'.  In general,
       the LSMR algorithm is supposed to be faster, however it has not
       been as thoroughly tested LSQR in this context.  Default is 'lsqr'

    path : str, optional
       The path where the HDF5 tables are stored.  Default is 'tables'

    minunc : float, optional
       The minimum value of the uncertainty.  This is to avoid divide by
       zero errors and infinite signal-to-noise.  Default is 1e-10
    """

    INT = np.uint64     # DO NOT CHANGE THIS
    FLOAT = np.float64  # DO NOT CHANGE THIS

    def __init__(self, extorders, mskorders=None, invmethod='lsqr', path='tables',
                 minunc=1e-10):

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
        return hasattr(self, 'A')

    def __str__(self):
        return f'Sparse matrix with {len(self)} elements'

    def build_matrix(self, grisms, sources, group=0):
        """
        Method to build a matrix

        Parameters
        ----------
        grisms : `WFSSCollection`
            The WFSS data to make into a matrix

        sources : `SourceCollection`
            The sources to build the matrix for

        group : int, optional
            The group ID.  Default is 0

        Notes
        -----
        1) This is the heart-and-soul of the multi-ended spectral extraction
        2) This operation can be very slow, but has been optimized to
           be efficient in its memory usage.

        """

        # do some quick checks
        if not self.extorders:
            LOGGER.warning('no valid extraction orders')
            return

        # we code it this way, so we can instantiate a matrix object
        # that contains all the defaults to update_headers.
        self.group = group

        # get some dimensionalities
        self.nfiles = len(grisms)  # grisms.nfiles
        # self.ndets=grisms.ndetectors
        self.nsources = len(sources)
        self.segids = list(sources.keys())

        # get the extraction parameters
        self.defpars = grisms.get_parameters()
        nwave = len(self.defpars)

        # compute the cumulative number of wavelength elements
        # this is to deal with uneven numbers of wavelength elements
        self.cumlam = [0]
        self.sedkeys = []         # to know which source/region goes with output
        for segid, source in sources.items():

            if hasattr(source, 'extpars'):
                nw = len(source.extpars)
            else:
                nw = nwave

            for sedkey, sedreg in source.items():
                self.cumlam.append(self.cumlam[-1]+nw)
                self.sedkeys.append(sedkey)
        self.cumlam = np.array(self.cumlam, dtype=int)

        # compute number of knowns (ie. total number of pixels) and
        # unknowns (ie. total number of wavelengths
        self.nknowns = grisms.npixels()     # total number of knowns
        self.nunknowns = self.cumlam[-1]    # total number of unknowns

        # record this value so we can make residuals
        self.imgdata = []

        # print a message
        LOGGER.info('Building a matrix\n' +
                    f'      Grism files: {self.nfiles}\n' +\
                    # f'      Grism dets:  {self.ndets}\n'+\
                    f'      Sources:     {self.nsources}\n' +\
                    f'      Unknowns:    {self.nunknowns}\n' +\
                    f'      Knowns:      {self.nknowns}')

        # lists to save the matrix content
        self.i = []     # knowns coordinate in matrix
        self.j = []     # unknowns coordinate in matrix
        self.aij = []   # matrix elements
        self.bi = []    # vector of known values

        # process each grism image
        self.detindex = 0
        for grism in tqdm(grisms, desc='Loading Grism Images', total=self.nfiles,
                          dynamic_ncols=True):

            # load the data for this grism image
            with PDTFile(grism, path=self.path, mode='r') as h5:

                for detname, detdata in grism.items():
                    h5.load_detector(detname)
                    flatfield = detdata.config.load_flatfield()

                    # record which pixels were used for this detector,
                    # this will be used later to get the unique pixels
                    self.xyg = []

                    # read the actual images
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

                    # make a good-pixel mask (gpx)
                    if detdata.config.bitmask:
                        gpx = np.bitwise_and(dqa, detdata.config.bitmask) == 0
                    else:
                        gpx = np.ones_like(dqa, dtype=bool)

                    # save the image dimensionalities
                    if not hasattr(self, 'imgdim'):
                        self.imgdim = tuple(reversed(sci.shape))

                    # updae msk for the mskorders
                    for ordname in self.mskorders:
                        h5.load_order(ordname)
                        for source in sources.values():
                            omt = h5.load_omt(source)
                            gpx[omt.get('y'), omt.get('x')] = False

                    # read each order to extract
                    for ordname in self.extorders:
                        h5.load_order(ordname)

                        # set a variable that will be updated in subroutine
                        self.srcindex = 0

                        # process each source
                        for source in sources.values():

                            for region in source:  # .spectralregions:
                                # read the table
                                tab = h5.load_odt(region)
                                if tab:
                                    self._parse_table(tab, source, unc, gpx,
                                                      ordname, detdata, flatfield)

                    # ok.  At this point, we've loaded all the data for
                    # a grism detector.  But let's verify that there are
                    # actually good pixels
                    if self.xyg:
                        # now must get the unique pixels from the observations
                        self.xyg = np.unique(self.xyg)

                        # compute the 2d pixel coordinates
                        xg, yg = np.unravel_index(self.xyg, detdata.naxis, order='F')

                        # grab the observations and weight by uncertainty
                        # (basically just the signal-to-noise)
                        bi = sci[yg, xg]/np.maximum(unc[yg, xg], self.minunc)

                        # save the results
                        self.bi.extend(bi)

                        # record the image info to make residual images later
                        self.imgdata.append((grism.dataset, detname))

                        # increment -- could just use len(self.imgdata)
                        self.detindex += 1

        # ok.  AT this point, we've loaded all the matrix data for the images
        #      but we'll need to do some juggling to best use the memory
        #      and computational resources

        # first we need to compress over the indices, which can be slow
        # so we'll print a message
        LOGGER.info('Compressing indices')
        self.icomp, self.iuniq = indices.compress(self.i)
        self.i.clear()          # delete the data

        self.jcomp, self.juniq = indices.compress(self.j)
        self.j.clear()

        # get some dimensionalities
        self.npix = len(self.iuniq)
        self.npar = len(self.juniq)
        dim = np.array([self.npix, self.npar], int)   # dimensionality of matrix

        # do a sanity check
        if np.amax(self.icomp) != len(self.bi)-1:
            msg = f'Matrix has invalid dimensionlity: ({self.npix}\u00d7{self.npar})'
            LOGGER.error(msg)
            raise RuntimeError(msg)

        # recast the vector
        self.bi = np.array(self.bi, dtype=float)

        # do some more checks
        if not self.overdetermined:
            LOGGER.warn(f"Underdetermined matrix: ({self.npix}\u00d7{self.npar}).\n" +
                        "Results will be dubious")

        # compute some extra things for the ragged arrays (ie. where
        # sources can have different number of spectral elements)
        # and there are tuples lurking to describe (segid,regid).  Will
        # call this as an extind for extraction index
        if self.nsources == 1:
            extind = np.zeros(self.npar, dtype=self.juniq.dtype)
        else:
            extind = np.digitize(self.juniq, self.cumlam)-1

        try:
            self.lamids = self.juniq-self.cumlam[extind]

        except BaseException:
            LOGGER.debug(self.npix, self.npar, self.nsources)
            raise RuntimeError("Matrix calulation has failed")

        # get the indices to do all reverse calculations
        self.ri = indices.reverse(extind)

        # compute the density of the matrix
        self.density = float(len(self.aij))/float(dim[0]*dim[1])

        # ok ok ok... NOW make a matrix
        self.A = coo_matrix((self.aij, (self.icomp, self.jcomp)), shape=dim)
        self.aij.clear()

        # compute the norm
        self.frob = linalg.norm(self.A)

        # make it a linear operator
        self.A = linalg.aslinearoperator(self.A)

        # compute the Frobenius norm and initialize the LCurve
        self.lcurve = LCurve(norm=self.frob)

        # now do a default damping target
        LOGGER.debug('damping target might be suspect for 2d objects')
        LOGGER.info("Building Damping Target")
        target = np.zeros(self.nunknowns, dtype=float)
        for segid, source in sources.items():

            for sedkey, region in source.items():
                extid = self.sedkeys.index(sedkey)

                g1 = self.ri[extid]
                g2 = self.lamids[g1]

                if hasattr(source, 'extpars'):
                    waves = source.extpars.wavelengths()
                else:
                    waves = self.defpars.wavelengths()

                target[g1] = region.sed(waves[g2])

        self.set_damping_target(target)

        # just update the printing...
        LOGGER.info('Finished building a matrix')

    def _parse_table(self, tab, source, unc, gpx, ordname, detdata, flatfield):
        """
        Helper function to parse data from a table

        Parameters
        ----------
        tab : `su.tables.HDF5Table`
            The table from which the data are taken

        source : `su.sources.Source`
            The source for which the table is to be read

        unc : `np.ndarray`
            The uncertainties from the image

        gpx : `np.ndarray`
            The good-pixel table

        ordname : str or int
            The name of the spectral order

        detdata : `DetectorConfig`
            The detector for which to parse the table

        flatfield : `su.core.wfss.config.FlatField`
            The flat field object to use


        Notes
        -----
        This is an internal helper function and should not really be
        explicitly called by a user.
        """

        # compute teh wavelengths
        wav = tab.wavelengths()
        x = tab.get('x', dtype=int)
        y = tab.get('y', dtype=int)
        val = tab.get('val')

        # remove any bad elements
        g = np.where(gpx[y, x])[0]

        x = x[g]
        y = y[g]
        val = val[g]

        # compute some things to rescale the matrix values
        # 1) sensitivity curve
        # 2) flat field
        # 3) pixel-area map (PAM)
        # 4) uncertainty
        sens = detdata.config[ordname].sensitivity(wav)*self.fluxscale
        flat = flatfield(x, y, wav)
        area = detdata.relative_pixelarea(x, y)
        uval = np.maximum(unc[y, x], self.minunc)     # set the uncertaint floor

        # rescale the tables' values
        val *= (sens*flat*area/uval)

        # group the wavelength indices according to the extraction settings
        # NB: the -1 is because digitize is 1-indexed, but arrays are
        #     0-indexed, so gotta shift it back
        if hasattr(source, 'extpars'):
            ll = np.digitize(wav, source.extpars.limits())-1
        else:
            ll = np.digitize(wav, self.defpars.limits())-1

        # compute the knowns index (called 'i' above).  First compute a
        # i-index offset based on number of pixels to add it here, but
        # will need to subtract it off below.
        ii0 = detdata.npixels()*self.detindex
        ii = np.ravel_multi_index((x, y), dims=detdata.naxis, order='F')+ii0

        # compute the unknowns index (called 'j' above), but gotta take
        # extra care here if there a different number of wavelength
        # elements per source.
        jj = ll+self.cumlam[self.srcindex]

        # decimate over the matrix coordinates --- meaning, sum over
        # repeated matrix indices (ii,jj).  This will leave only the
        # unique pairs of (i,j), hence the 'u' designation
        vu, iu, ju = indices.decimate(val, ii, jj, dims=(self.nknowns, self.nunknowns))

        # get the grism-pixel coordinates (subtracting the offset as promised)
        xygu = indices.uniq(iu-ii0)

        # save all these values and increment a global counter
        self.i.extend(list(iu))
        self.j.extend(list(ju))
        self.aij.extend(list(vu))
        self.xyg.extend(list(xygu))
        self.srcindex += 1

    def compute_model(self, xj):
        """
        Method to compute the model, which is given as a matrix-product
        between this object and a vector of spectral values

        Parameters
        ----------
        xj : `np.ndarray`
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
        if self:

            bi = self.A.matvec(xj)

            npix = self.imgdim[0]*self.imgdim[1]

            imgindices, pixindices = np.divmod(self.iuniq, npix)

            dtype = [('x', np.uint16), ('y', np.uint16), ('v', np.float64)]

            out = {image: dict() for (image, detname) in self.imgdata}
            for imgindex, (image, detname) in enumerate(self.imgdata):

                g = np.where(imgindices == imgindex)

                y, x = np.divmod(pixindices[g], self.imgdim[0])

                out[image][detname] = np.array(list(zip(x, y, bi[g])), dtype=dtype)

        else:
            out = None

        return out

    def set_damping_target(self, target):
        """
        Method to set the damping target

        Parameters
        ----------
        target : `np.ndarray`
           The damping target as a concatenated array of spectra

        """
        if self:
            self.target = target/self.fluxscale
            self.bi -= self.A.matvec(self.target)
        else:
            LOGGER.warning("Cannot set damping target until matrix is built")

    # some flagging parameters

    @property
    def underdetermined(self):
        return self.npar > self.npix

    @property
    def determined(self):
        return self.npar == self.npix

    @property
    def overdetermined(self):
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

    # stuff for inversions
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

        r = linalg.lsqr(self.A, self.bi, damp=damp, calc_var=True, **kwargs)
        (x, istop, itn, r1norm, r2norm, anorm, acond, arnorm, xnorm, var) = r

        r = Result('LSQR', x, istop, itn, r1norm, r2norm, anorm, acond, arnorm, xnorm,
                   np.sqrt(var), damp)

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

        (x, istop, itn, norm, arnorm, anorm, acond, xnorm) = linalg.lsmr(self.A, self.bi,
                                                                         damp=damp, **kwargs)
        r1norm = arnorm
        r2norm = np.sqrt(arnorm**2+(damp*xnorm)**2)
        unc = np.full_like(x, np.nan)
        r = Result('LSMR', x, istop, itn, r1norm, r2norm, anorm, acond, arnorm, xnorm,
                   unc, damp)
        return r

    def compute_uncertainty(self):
        """
        Function to compute the uncertainties

        Notes
        -----
        The presence of damping means that the inversion algorithms
        (either LSQR or LSMR) will generate uncertainties that are too small
        as they are estimating diagonals of an altered matrix (Abar).
        See the documentation for LSQR and the `calc_var` parameter
        for a description of that altered matrix.  Therefore, this method
        estimates the diagonal elements of AT.A and inverting them. This is
        not explicitly correct, as we need the diagonal elements of (AT.A)-1.

        """
        LOGGER.warning("Uncertainties are not exact in `matrix.py`")

        aij2, j = indices.decimate(self.A.A.data*self.A.A.data,
                                   self.juniq[self.jcomp])

        # these two lines to avoid divide-by-zero errors
        unc = np.full_like(aij2, np.inf)
        np.divide(1., np.sqrt(aij2), out=unc, where=(aij2 != 0))

        return unc

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
            LOGGER.warning("Cannot update header, matrix unitialized")
            return

        hdr['NNZ'] = (len(self), 'number of non-zero matrix elements')

        hdr['NPIX'] = (self.npix, 'number of pixels analyzed (ie. knowns)')
        hdr['NPAR'] = (self.npar, 'number of parameters measured (ie. unknowns)')
        hdr['NDOF'] = (self.npix-self.npar, 'number of degrees of freedom')
        hdr['DENSITY'] = (self.density, 'fraction of non-zero elements')
        hdr['FROBNORM'] = (self.frob, 'Frobenius norm')

        hdr['INVMETH'] = (self.invmethod, 'Method to compute inverse')
        hdr['MINUNC'] = (self.minunc, 'minimum uncertainty in e/s')
        hdr['EXTORDER'] = (','.join(self.extorders), 'extraction orders')
        hdr['MSKORDER'] = (','.join(self.mskorders), 'masked orders')
        headers.add_stanza(hdr, 'Matrix Properties', before='NNZ')

    # --------------------------------------------------------------------
    #
    #                       helper functions
    #
    # --------------------------------------------------------------------
    # @staticmethod
    # def tuplize(d):
    #    if d is not None:
    #        if not isinstance(d,(tuple,list)):
    #            d=(d,)
    #    else:
    #        d=()
    #    return d

    # def retype(self,k,dtype):
    #    if dtype is np.nan:
    #        if hasattr(self,k):
    #            return getattr(self,k)
    #        else:
    #            return np.nan
    #    else:
    #        if hasattr(self,k):
    #            return dtype(getattr(self,k))
    #        else:
    #            return dtype(0)
