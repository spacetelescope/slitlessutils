import os

import numpy as np
from astropy.io import fits
from astropy.stats import SigmaClip

from .....config import SUFFIXES, Config
from .....info import __code__
from .....logger import LOGGER
from ....tables import PDTFile
from ....utilities import headers, indices  # , as_iterable
from ...module import Module
from .contamination import Contamination


class Single(Module):
    """
    Module to implement single-ended extraction, which is conceptually
    similar to the Horne (1986) approach.

    inherits from `su.core.modules.Module`

    Parameters
    ----------
    extorders : str, or iterable
       The name of the spectral orders to subtract.  This is likely a scalar
       string, but could be an iterable of strings to extract multiple
       spectral orders.

    savecont : bool, optional
       A flag to save the contamination.  Default is False

    nsigmaclip : int, float, or 2-element iterable, optional
       The number of sigma to clip bad pixels.  If this is a scalar value,
       then it represents the low/high thresholds.  If this is a 2-element
       iterable, then the two values are for independent low/high
       (respectively) values.  See also `astropy.stats.SigmaClip()`
       Default is 3.

    maxiters : int, optional
       The maximum number in the sigma clipping.  See also
       `astropy.stats.SigmaClip()` Default is 10.

    root : str, optional
        The output fits file root name, which is the primary output of 1d
        spectra.  If this is set to None, then "slitlessutils" is used and
        hence will be overwritten following subsequent runs, so it is advised
        to set this.  Default is None.

    outpath : str, optional
        A path where output products will be written.  If set to None or invalid,
        then CWD will be used.  If set to valid str, then path will be created
        (if possible), then used.  If still not valid, then will use CWD.
        Default is None.

    writecsv : bool, optional
        Flag to write light-weight boolean files.  Default is True.

    kwargs : dict, optional
       Dictionary of additional arguments passed to `su.core.modules.Module`

    Notes
    -----
    The canonical Horne (1986) approach is very tricky to implement for
    WFSS data, as there are many sums over the cross dispersion profile.
    But many WFSS instruments have a cross-dispersion profile that is
    not parallel to either the x- or y-axes, or even has many twists
    and turns.  Therefore those sums can be tedious to define, let alone
    compute.  This approach is conceptually similar, but framed differently.

    Example
    -------
    The main entry point is through the `__call__()` method in the `Module`
    parent class.

    >>> import slitlessutils as su  # doctest: +SKIP
    >>> ext = su.modules.Single('+1', mskorders=None, root=ROOT)  # doctest: +SKIP
    >>> res = ext(data,sources)  # doctest: +SKIP

    """

    DESCRIPTION = 'Extracting (single)'

    DTYPE = [('lamb', np.float32),
             ('flam', np.float32),
             ('func', np.float32),
             ('cont', np.float32),
             ('npix', np.uint16)]

    FILETYPE = '1d spectra'

    def __init__(self, extorders, mskorders='all', savecont=False, root=None,
                 outpath=None, writecsv=True, nsigmaclip=3, maxiters=10,
                 **kwargs):
        Module.__init__(self, self.extract, postfunc=self._combine, **kwargs)

        # self.extorders=as_iterable(extorders)
        # specify the orders
        self.extorders = self.as_set(extorders)
        self.mskorders = self.as_set(mskorders)
        if self.extorders.intersection(self.mskorders):
            msg = 'Masking an extraction order is not ok'
            LOGGER.warning(msg)
            raise RuntimeError(msg)

        # set up for the contamination
        self.contamination = Contamination(mskorders)
        self.savecont = savecont

        # for sigma-clipping of the data
        if isinstance(nsigmaclip, (int, float)):
            self.sigclip = SigmaClip(sigma=nsigmaclip, maxiters=maxiters)

        # set output file path
        self.outpath = self.set_outpath(outpath)
        self.writecsv = writecsv
        if self.writecsv:
            self.csvpath = self.outpath + 'CSVFILES' + os.sep
            if not os.path.exists(self.csvpath):
                try:
                    os.mkdir(self.csvpath)
                except FileNotFoundError:
                    LOGGER.warning(f'Cannot make CSVPATH {self.csvpath},'
                                   ' so no CSV files will be written')
                    self.writecsv = False
                    self.csvpath = None

        # output file names
        if root is None:
            root = __code__
        self.filename = os.path.join(self.outpath, f"{root}_{SUFFIXES[self.FILETYPE]}.fits")

        if self.savecont:
            msg = "no support to write contamination images yet"
            LOGGER.warning(msg)
            # self.savecont = False

    @staticmethod
    def set_outpath(path):
        """
        Static method to implement rules for setting output path.

        The rules are:
        1) if path is not a valid string, then use CWD
        2) if path is valid string then
            a) if path is an existing path
                i) if path is writable, then return that
                ii) if path is not writable, then return CWD
            b) if path is not an existing path try making the dir.
                i) if can make, then return that
                ii) if cannot make, then return CWD

        Parameters
        ----------
        path : str
            Path to consider as the output path

        Returns
        -------
        outpath : str
            A qualified, valid path
        """

        outpath = os.getcwd()
        if isinstance(path, str):
            if os.path.isdir(path):
                if os.access(path, os.W_OK):
                    # if here, then has valid data type and a writable path
                    outpath = path
            else:
                # passed a str, but not a valid path.  try making the dir?
                try:
                    os.mkdir(path)
                    outpath = path
                except FileNotFoundError:
                    pass
        if outpath[-1] != os.sep:
            outpath += os.sep

        return outpath

    def _combine(self, results, data, sources, **kwargs):
        """
        Method to aggregate the spectra for each source from the separate
        WFSS exposures.

        Parameters
        ----------
        results : dict
            A dictionary of dictionaries, one for each segid.

        data : `su.core.wfss.WFSSCollection`
            The collection of WFSS images

        sources : `su.core.sources.SourceColection`
            The collection of sources

        Notes
        -----
        This is likely not to be directly called.
        """

        LOGGER.info("Combining the 1d spectra")

        # suffix of all output data
        prekwargs = {'filetype': (self.FILETYPE, 'contents of this file'),
                     'exttype': ('single', 'method for extraction'),
                     'savecont': (self.savecont, 'Was contamination saved?')}

        # grab the default extraction parameters
        defpars = data.get_parameters()

        # the below may not work if a source isn't in the first result catalog?

        # aggregate the results
        result = results.pop(0)
        segids = list(result.keys())
        for segid in segids:
            result[segid]['file'] = [0]*len(result[segid]['flam'])

        i = 1
        while results:
            r = results.pop(0)
            for segid, res in r.items():
                result[segid]['flam'].extend(res['flam'])
                result[segid]['func'].extend(res['func'])
                result[segid]['wave'].extend(res['wave'])
                result[segid]['cont'].extend(res['cont'])
                result[segid]['file'].extend([i]*len(res['flam']))
                i += 1

        # get all the objects in this output catalog
        segids = list(result.keys())

        # get the flux funits
        config = Config()

        # make an output structure
        hdul = fits.HDUList()

        # make a primary header and fill it with useful info
        phdu = fits.PrimaryHDU()
        headers.add_preamble(phdu.header, **prekwargs)  # basic preamble
        headers.add_software_log(phdu.header)           # about software
        data.update_header(phdu.header)                 # about WFSS images
        sources.update_header(phdu.header)              # about the sources
        config.update_header(phdu.header)               # global config info
        self.contamination.update_header(phdu.header)   # about the contamination

        # if there is a sigma clipper
        phdu.header.set('SIGCLIP', value=hasattr(self, 'sigclip'),
                        comment='Sigma-clip weighting?')
        if phdu.header['SIGCLIP']:
            phdu.header.set('NSIGMA', value=self.sigclip.sigma,
                            comment='number of sigma for clipping')
            phdu.header.set('MAXITER', value=self.sigclip.maxiters,
                            comment='maximum number of iterations')
        headers.add_stanza(phdu.header, 'Spectral Combination', before='SIGCLIP')
        super().update_header(phdu.header)            # about the CPU settings

        # add to the HDUList
        hdul.append(phdu)

        # now process each object
        for segid in segids:
            res = result.pop(segid)

            source = sources[segid]

            # set the spectral extraction parameters
            pars = defpars.update_pars(source[0])

            # convert to np.array to do calculation
            flam = np.array(res['flam'])
            func = np.array(res['func'])
            wave = np.array(res['wave'])
            cont = np.array(res['cont'])
            file = np.array(res['file'])

            # find the spectral bins that these elements belong to
            lamb = pars.indices(wave)

            # number of wavelength elements
            nwave = len(pars)
            g = np.where((0 < lamb) & (lamb < nwave))[0]

            # get good values
            flam = flam[g]        # spectrum
            func = func[g]        # uncertainties
            cont = cont[g]        # contamination model
            wave = wave[g]        # wavelengths in A
            lamb = lamb[g]        # wavelength indices
            file = file[g]

            # get reverse elements
            ri = indices.reverse(lamb)

            # make an output data structure
            out = np.zeros(nwave, dtype=self.DTYPE)
            out['lamb'] = pars.wavelengths()           # PRISM CHANGE
            out['flam'] = np.nan
            out['func'] = np.nan
            out['cont'] = np.nan
            out['npix'] = 0

            # compute weighted-averages over the bins
            for ind, g in ri.items():
                val = flam[g]
                unc = func[g]
                wht = 1./unc**2      # initialize weights as inverse variance
                cnt = cont[g]

                # update the weights with sigma clipping if requested
                # if hasattr(self, 'sigclip'):
                #    masked = self.sigclip(val, masked=True, copy=True)
                #    wht *= (1-np.ma.getmask(masked).astype(float))

                # if a pixel has an infinite cont, then ignore it
                # wht=np.where(np.isinf(cnt),0.,wht)

                # update weights for the contamination
                # gc=np.where(cnt>0)[0]
                # if len(gc)>0:
                #    frac=cnt[gc]/val[gc]
                #    wht[gc]/=frac**2

                # weighted STDEV as the error?
                # err = np.sqrt(np.average((val-ave)**2, weights=wht))

                # compute sum of weights (used later)
                nwht = np.count_nonzero(wht)
                if nwht > 0:
                    wsum = np.sum(wht)

                    # iteratively clip points
                    # for it in range(1):
                    #    wsum = np.sum(wht)
                    #    ave = np.sum(wht*val)/wsum
                    #    delta = np.abs(val - ave)/unc
                    #    bad = np.where(np.abs(delta)>3.)[0]
                    #    wht[bad] = 0.0

                    # save the results
                    out['flam'][ind] = np.sum(wht*val)/wsum
                    out['func'][ind] = 1./np.sqrt(wsum)
                    out['cont'][ind] = np.sum(wht*cnt)/wsum
                    out['npix'][ind] = nwht

            # update the source
            source[0].sed.reset(out['lamb'], out['flam'], func=out['func'],
                                cont=out['cont'], npix=out['npix'])

            # write out 1d SEDs?
            if self.writecsv:
                source.write_seds(filetype='csv', path=self.csvpath)

            # store source in the OUTPUT file
            hdul.append(source.as_HDU())

        # get basename of output file
        LOGGER.info(f'Writing: {self.filename}')
        hdul.writeto(self.filename, overwrite=True)

    @staticmethod
    def apply_bitmask(dqa, *args, bitmask=None):
        """
        Staticmethod to remove elements which have bad values in the data
        quality arrays (DQAs).

        Parameters
        ----------
        dqa : `np.ndarray`
            the data quality array as a np array

        *args : tuple of `np.ndarray`s
            the data to remove the bad pixels from

        bitmask : int or None
            The bitmask to consider bad pixels.  If None, then no flagging is done.
            Default is None.

        Returns
        -------
        outs : a tuple of lists of the same datatype as the *args
        """

        if bitmask:
            g = np.where(np.bitwise_and(dqa, bitmask) == 0)[0]
            if g.size > 0:
                out = tuple(a[g] for a in args)
            else:
                LOGGER.warning(f'Bitmask ({bitmask}) removes all pixels')
                out = tuple([] for a in args)
            if len(args) == 1:
                out = out[0]
            return out
        else:
            return args

    def extract(self, data, sources, **kwargs):
        """
        Method to do the single-ended spectral extraction

        Parameters
        ----------
        data : `su.core.wfss.WFSS`
            A WFSS image to extract

        sources : `su.core.sources.SourceColection`
            The collection of sources

        cartesian : bool, optional
            Flag that extraction should be exactly along pixel rows/columns.
            This is the behavior of aXe.  Default is True

        profile : str, optional
            Reset the cross dispersion weights.  Can be:

            'forward': use the forward model
            'uniform': use a uniform/flat cross dispersion profile, which is
                effectively just a box-extraction
            'data': use the WFSS data as the profile, which is
                effectively the Horne 1986 setting without the wavelength-
                dependent smoothing

            Default is 'forward'

        Returns
        -------
        results : dict
            This is a dictionary of dictionaries that get passed to
            `self._combine`

        Notes
        -----
        This is likely not to be directly called.

        """

        # scale of the output fluxes
        fluxscale = Config().fluxscale

        # create a data structure to hold the results
        results = {segid: {'wave': [], 'flam': [], 'func': [], 'cont': [], 'dwav': []} for segid in
                   sources.keys()}

        # sort out optional inputs
        cartesian = kwargs.get('cartesian', True)
        profile = kwargs.get('profile', 'data').lower()

        # padding for the contamination models
        padx = (5, 5)
        pady = (5, 5)

        # contamination image
        if self.savecont:
            chdul = fits.HDUList()
            phdr = fits.PrimaryHDU()

            LOGGER.debug("add stuff to contam images' primary headers")
            chdul.append(phdr)

        # grab the default extraction parameters
        defpars = data.get_parameters()

        # open the PDT for reading
        with PDTFile(data, path=self.path, mode='r') as h5:

            for detname, detdata in data.items():
                h5.load_detector(detname)
                bitmask = detdata.config.bitmask

                # load the data
                phdr = detdata.primaryheader()
                sci, hdr = detdata.readfits('science', header=True)
                unc = detdata.readfits('uncertainty')
                dqa = detdata.readfits('dataquality')

                # get the flat field
                flatfield = detdata.config.load_flatfield(**kwargs)

                # convert units to a rate
                if hdr['BUNIT'].lower() in ('electron', 'electrons', 'e', 'e-'):
                    time = phdr['EXPTIME']
                    sci /= time
                    unc /= time

                # initialize the contamination modeling
                self.contamination.make_model(sources, h5)

                # process each order
                orders = self.extorders if self.extorders else detdata.orders
                for ordname in orders:
                    h5.load_order(ordname)

                    # process each source
                    for segid, source in sources.items():

                        odt = h5.load_odt(source)
                        if odt:
                            # get contents of the table
                            xg = odt.get('x')
                            yg = odt.get('y')
                            val = odt.get('val')
                            # lam = odt.get('lam')
                            wav = odt.wavelengths()
                            dwav = np.zeros_like(val, dtype=float)

                            # set the spectral extraction parameters
                            pars = defpars.update_pars(source[0])
                            wavelengths = pars.wavelengths()

                            # limits = pars.limits()
                            # lam = pars.indices(wav)

                            # get minimum bounding box for lops
                            x0, x1, y0, y1 = odt.bounding_box()

                            # if requested, compute the contamination model:
                            if self.contamination:
                                # grow bounding boxes for contamination models
                                bbx = (np.maximum(x0-padx[0], 0),
                                       np.minimum(x1+padx[1], detdata.naxis[0]-1))
                                bby = (np.maximum(y0-pady[0], 0),
                                       np.minimum(y1+pady[1], detdata.naxis[1]-1))

                                # get the contam data
                                chdu = self.contamination(segid, ordname, h5, sources,
                                                          detdata, flatfield, bbx=bbx, bby=bby)

                                if self.savecont:
                                    chdul.append(chdu)

                                # get the offsets in the contamination image
                                xoff = -chdu.header.get('LTV1', 0)
                                yoff = -chdu.header.get('LTV2', 0)

                            if cartesian:
                                # for making the residuals
                                # mod = np.zeros_like(sci)
                                for x in range(x0, x1+1, 1):
                                    g = np.where(x == xg)[0]
                                    xx = xg[g]
                                    yy = yg[g]
                                    ww = wav[g]
                                    vv = val[g]
                                    dw = dwav[g]

                                    # apply the bitmask
                                    xx, yy, vv, ww, dw = self.apply_bitmask(
                                        dqa[yy, xx], xx, yy, vv, ww, dw, bitmask=bitmask)

                                    ww, _y = indices.decimate(ww*vv, yy)
                                    vv, yy = indices.decimate(vv, yy)
                                    ww /= vv
                                    xx = np.full_like(yy, x, dtype=int)

                                    # the average wavelength in this column
                                    wave = np.average(ww, weights=vv)

                                    # compute the calibrations
                                    disp = detdata.config[ordname].dispersion(
                                        *source.xyc, wavelength=ww)
                                    sens = detdata.config[ordname].sensitivity(ww)
                                    flat = flatfield(xx, yy, ww)
                                    area = detdata.relative_pixelarea(xx, yy)
                                    den = flat*area*sens*fluxscale*disp

                                    # apply calibrations
                                    ss = sci[yy, xx]/den
                                    uu = unc[yy, xx]/den

                                    # set the cross dispersion profile
                                    if profile == 'uniform':
                                        # simple summing over pixels
                                        flam = np.sum(ss)
                                        func = np.sqrt(np.sum(uu**2))
                                    else:
                                        # weighting by profile (a la Horne)
                                        if profile == 'forward':
                                            prof = vv.copy()
                                        elif profile == 'data':
                                            prof = np.maximum(sci[yy, x], 0.)
                                        else:
                                            msg = f'Profile setting ({profile}) is invalid'
                                            LOGGER.error(msg)
                                            raise RuntimeError(msg)

                                        # normalize the profile
                                        prof /= np.sum(prof)

                                        wht = prof/uu**2
                                        norm = np.sum(prof*wht)
                                        flam = np.sum(ss*wht)/norm
                                        func = np.sqrt(np.sum(prof)/norm)

                                    # this can happen if sens==0
                                    if np.isnan(flam):
                                        flam = 0.0

                                    # update the model
                                    # mod[yy, x] = flam*den*prof

                                    # compute contamination
                                    if self.contamination:
                                        cc = chdu.data[yy-yoff, xx-xoff]/den
                                        if profile == 'uniform':
                                            cont = np.sum(cc)
                                        else:
                                            cont = np.sum(cc*wht)/norm
                                    else:
                                        cont = 0.0

                                    # save the results
                                    results[segid]['flam'].append(flam)
                                    results[segid]['func'].append(func)
                                    results[segid]['wave'].append(wave)
                                    results[segid]['dwav'].append(0.0)
                                    results[segid]['cont'].append(cont)

                            else:

                                di = 0.5
                                for ind, wave in enumerate(wavelengths):
                                    w0, w1 = pars(ind-di), pars(ind+di)
                                    g = np.where((w0 <= wav) & (wav < w1))[0]
                                    # for ind, g in ri.items():
                                    xx = xg[g]
                                    yy = yg[g]
                                    ww = wav[g]
                                    vv = val[g]
                                    dw = dwav[g]

                                    # apply the bitmask
                                    xx, yy, vv, ww, dw = self.apply_bitmask(
                                        dqa[yy, xx], xx, yy, vv, ww, dw, bitmask=bitmask)

                                    # compute cross dispersion direction
                                    # t = detdata.config[ordname].displ.invert(source.xyc[0],
                                    # source.xyc[1], thiswave)
                                    # dxdt = detdata.config[ordname].dispx.deriv(source.xyc[0],
                                    # source.xyc[1], t)
                                    # dydt = detdata.config[ordname].dispy.deriv(source.xyc[0],
                                    # source.xyc[1], t)
                                    # drdt = np.sqrt(dxdt**2 + dydt**2)

                                    # dispx = dxdt/drdt
                                    # dispy = dydt/drdt

                                    # cross dispersion vector (orthogonal, by construction)
                                    # crossx = -dy
                                    # crossy = +dx

                                    # center of region
                                    # x0 = np.average(xx, weights=vv)
                                    # y0 = np.average(yy, weights=vv)

                                    # distances to the cross dispersion vector
                                    # D = np.abs(crossx*(yy-y0)-crossy*(xx-x0))
                                    # g = np.where(D<2)[0]
                                    # if g.size>0:
                                    #     xx = xx[g]
                                    #     yy = yy[g]
                                    #     vv = vv[g]
                                    #     ww = ww[g]
                                    #     dw = dw[g]

                                    # measure size
                                    # X = np.average(xx, weights=vv)
                                    # Y = np.average(yy, weights=vv)
                                    # X2 = np.average((xx-X)**2, weights=vv)
                                    # Y2 = np.average((yy-Y)**2, weights=vv)
                                    # XY = np.average((xx-X)*(yy-Y), weights=vv)
                                    # DIF = (X2-Y2)/2.
                                    # AVE = (X2+Y2)/2.
                                    # DET = DIF**2 + XY**2
                                    # A = np.sqrt(AVE + np.sqrt(DET))
                                    # B = np.sqrt(AVE - np.sqrt(DET))
                                    # SIG = np.sqrt(A*B)
                                    # FWHM = SIG*2.35

                                    # decimate over unique pixels
                                    ww, _x, _y = indices.decimate(ww*vv, xx, yy)
                                    vv, xx, yy = indices.decimate(vv, xx, yy)
                                    ww /= vv

                                    disp = detdata.config[ordname].dispersion(
                                        *source.xyc, wavelength=wave)
                                    sens = detdata.config[ordname].sensitivity(ww)
                                    flat = flatfield(xx, yy, ww)
                                    area = detdata.relative_pixelarea(xx, yy)
                                    den = flat*area*sens*fluxscale*disp
                                    ss = sci[yy, xx]/den
                                    uu = unc[yy, xx]/den

                                    # using forward model weights
                                    prof = np.maximum(vv, 0.)
                                    prof = prof/np.sum(prof)
                                    wht = prof/uu**2
                                    norm = np.sum(prof*wht)
                                    flam = np.sum(ss*wht)/norm
                                    func = np.sqrt(np.sum(prof)/norm)

                                    # if doing contamination
                                    if self.contamination:
                                        cc = chdu.data[yy-yoff, xx-xoff]
                                        cont = np.sum(cc*wht)/norm
                                    else:
                                        cont = 0.0

                                    # store the results
                                    results[segid]['flam'].append(flam)
                                    results[segid]['func'].append(func)
                                    results[segid]['wave'].append(wave)
                                    results[segid]['dwav'].append(0.0)
                                    results[segid]['cont'].append(cont)

        if self.savecont:
            chdul.writeto(f'{data.dataset}_cont.fits', overwrite=True)

        return results
