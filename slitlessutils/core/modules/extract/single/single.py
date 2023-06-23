from astropy.io import fits
from astropy.stats import SigmaClip
import numpy as np

from .....config import Config, SUFFIXES
from .contamination import Contamination
from .....info import __code__
from .....logger import LOGGER
from ...module import Module
from ....tables import PDTFile
from ....utilities import indices, headers  # , as_iterable


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

    >>> import slitlessutils as su
    >>> ext = su.modules.Single('+1',mskorders=None,root=ROOT)
    >>> res = ext(data,sources)

    """

    DESCRIPTION = 'Extracting (single)'

    DTYPE = [('lamb', np.float32),
             ('flam', np.float32),
             ('func', np.float32),
             ('cont', np.float32),
             ('npix', np.uint16)]

    FILETYPE = '1d spectra'

    def __init__(self, extorders, mskorders='all', savecont=False,
                 nsigmaclip=3, maxiters=10, **kwargs):
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

        # output file names
        root = kwargs.get('root', __code__)
        self.filename = f"{root}_{SUFFIXES[self.FILETYPE]}.fits"

        if self.savecont:
            raise NotImplementedError("no support to write contam images yet")

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

        # grab the default extraction parmaeters
        defpars = data.get_parameters()

        # the below may not work if a source isn't in the first result catalog?

        # aggregate the results
        result = results.pop(0)
        while results:
            r = results.pop(0)
            for segid, res in r.items():
                result[segid]['flam'].extend(res['flam'])
                result[segid]['func'].extend(res['func'])
                result[segid]['wave'].extend(res['wave'])
                result[segid]['cont'].extend(res['cont'])

        # get all the objects in this output catalog
        segids = list(result.keys())

        # get the flux funits
        config = Config()
        # funits = f'{config.fluxscale} {config.fluxunits}'

        # make an output structure
        hdul = fits.HDUList()

        # make a primary header and fill it with useful info
        phdu = fits.PrimaryHDU()
        headers.add_preamble(phdu.header, **prekwargs)  # basic preamble
        headers.add_software_log(phdu.header)         # about software
        data.update_header(phdu.header)           # about WFSS images
        sources.update_header(phdu.header)            # about the sources
        config.update_header(phdu.header)             # global config info
        self.contamination.update_header(phdu.header)  # about the contamination

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

            # get reverse elements
            ri = indices.reverse(lamb)

            # make an output data structure
            out = np.zeros(nwave, dtype=self.DTYPE)
            out['lamb'] = pars.wavelengths()
            out['flam'] = np.nan
            out['func'] = np.nan
            out['cont'] = np.nan
            out['npix'] = 0

            # compute weighted-averages over the bins
            for ind, g in ri.items():
                val = flam[g]
                wht = 1./func[g]**2
                cnt = cont[g]

                # update the weights with sigma clipping if requested
                if hasattr(self, 'sigclip'):
                    masked = self.sigclip(val, masked=True, copy=True)
                    wht *= (1-np.ma.getmask(masked).astype(float))

                # if a pixel has an infinite cont, then ignore it
                # wht=np.where(np.isinf(cnt),0.,wht)

                # update weights for the contamination
                # gc=np.where(cnt>0)[0]
                # if len(gc)>0:
                #    frac=cnt[gc]/val[gc]
                #    wht[gc]/=frac**2

                # ave = np.average(val, weights=wht)
                # err = np.sqrt(np.average((val-ave)**2, weights=wht))

                # compute sum of weights (used later)
                wsum = np.sum(wht)

                # save the results
                out['flam'][ind] = np.sum(wht*val)/wsum
                out['func'][ind] = 1./np.sqrt(np.sum(wht))
                out['cont'][ind] = np.sum(wht*cnt)/wsum
                out['npix'][ind] = np.count_nonzero(wht)

            # update the source
            source[0].sed.reset(out['lamb'], out['flam'], func=out['func'],
                                cont=out['cont'], npix=out['npix'])

            # write out 1d SEDs
            source.write_seds(filetype='csv')

            # store source in the OUTPUT file
            hdul.append(source.as_HDU())

        # get basename of output file
        LOGGER.info(f'Writing: {self.filename}')
        hdul.writeto(self.filename, overwrite=True)

    def extract(self, data, sources, **kwargs):
        """
        Method to do the single-ended spectral extraction

        Parameters
        ----------
        data : `su.core.wfss.WFSSCollection`
            The collection of WFSS images

        sources : `su.core.sources.SourceColection`
            The collection of sources

        Returns
        -------
        results : dict
            This is a dictionary of diciontaries that get passed to
            `self._combine`

        Notes
        -----
        This is likely not to be directly called.

        """

        # scale of the output fluxes
        fluxscale = Config().fluxscale

        # create a data structure to hold the results
        results = {segid: {'wave': [], 'flam': [], 'func': [], 'cont': []} for segid in
                   sources.keys()}

        # contamination image
        # if self.savecont:
        #    hdul=fits.HDUList()
        #    phdr=fits.PrimaryHDU()
        #    LOGGER.debug("add stuff to contam images' primary headers")
        #    hdul.append(phdr)

        hdul = fits.HDUList()

        # open the PDT for reading
        with PDTFile(data, path=self.path, mode='r') as h5:

            for detname, detdata in data.items():
                h5.load_detector(detname)
                dims = detdata.naxis

                # load the data
                phdr = detdata.primaryheader()
                sci, hdr = detdata.readfits('science', header=True)
                unc = detdata.readfits('uncertainty')
                # dqa = detdata.readfits('dataquality')

                # get the flat field
                flatfield = detdata.config.load_flatfield(**kwargs)

                # convert units to a rate
                if hdr['BUNIT'].lower() in ('electron', 'electrons', 'e', 'e-'):
                    time = phdr['EXPTIME']
                    sci /= time
                    unc /= time

                # initialize the contamination modeling
                self.contamination.make_model(sources, h5)
                # cont=Contam1(sources,h5)
                # contmodel=self.contamination.initialize_model(hdr)
                # contmodel.update_model_from_hdf5(sources,h5)

                # process each order
                orders = self.extorders if self.extorders else detdata.orders
                for ordname in orders:
                    h5.load_order(ordname)

                    # process each source
                    for segid, source in sources.items():

                        # load the object-profile table
                        opt = h5.load_opt(source)
                        if opt:

                            # get contents of the table
                            xg = opt.get('x')
                            yg = opt.get('y')
                            val = opt.get('val')
                            wav = opt.get('wav')

                            # get the contamination model
                            # cnt=contmodel.weight(xg,yg)

                            # find the unique grism-image pixels and compute
                            # the average wavelength per pixel
                            vv, xx, yy = indices.decimate(val, xg, yg, dims=dims)
                            ww, _x, _y = indices.decimate(wav*val, xg, yg, dims=dims)

                            # compute the average wavelngth in the pixel
                            ww /= vv

                            # get data quality
                            # d = dqa[yy, xx]

                            # TODO: GOTTA REMOVE PIXELS FROM DQA

                            # compute correction factors
                            sens = detdata.config[ordname].sensitivity(ww)
                            flat = flatfield(xx, yy, ww)
                            area = detdata.relative_pixelarea(xx, yy)

                            # aggregate factors
                            den = vv*sens*flat*area*fluxscale
                            g = np.where((den > 0))[0]

                            # only keep terms that are strictly positive
                            if g.size > 0:
                                xx = xx[g]
                                yy = yy[g]
                                vv = vv[g]
                                ww = ww[g]
                                den = den[g]

                                # grab the science/uncertainty
                                s = sci[yy, xx]
                                u = unc[yy, xx]

                                # null out the contamination
                                c = np.zeros_like(s, dtype=float)

                                # if doing contamination
                                if self.contamination:

                                    # get the contamination
                                    hdu = self.contamination(segid, ordname, h5,
                                                             sources, detdata,
                                                             flatfield)
                                    # write the contamination to a file?
                                    if self.savecont:
                                        hdul.append(hdu)

                                    # grab dimensionalities of the cont image
                                    x0 = -hdu.header['LTV1']
                                    nx = hdu.header['NAXIS1']

                                    y0 = -hdu.header['LTV2']
                                    ny = hdu.header['NAXIS2']

                                    # only work with pixels in this subimage
                                    gg = np.where((x0 < xx) & (xx <= x0+nx) &
                                                  (y0 < yy) & (yy <= y0+ny))[0]
                                    if gg.size > 0:
                                        # save the contamination
                                        c[gg] = hdu.data[yy[gg]-y0, xx[gg]-x0]

                                # store the results
                                results[segid]['flam'].extend(s/den)
                                results[segid]['func'].extend(u/den)
                                results[segid]['wave'].extend(ww)
                                results[segid]['cont'].extend(c/den)

        if self.savecont:
            hdul.writeto(f'{data.dataset}_cont.fits', overwrite=True)

        return results
