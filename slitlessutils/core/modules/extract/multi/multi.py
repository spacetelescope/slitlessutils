import getpass

from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
from tqdm import tqdm

from .....config import SUFFIXES, Config
from .....logger import LOGGER
from ....utilities import as_iterable, get_metadata, headers, compression
from ...group import GroupCollection
from ...module import Module
from .matrix import Matrix
from .optimizer import optimizer


class Multi(Module):
    """
    Module to implement multi-ended spectral extraction, as detailed in
    Ryan, Casertano, & Pirzkal (2018)

    inherits from `su.core.modules.Module`

    Parameters
    ----------
    extorders : str, or iterable
       The name of the spectral orders to subtract.  This is likely a scalar
       string, but could be an iterable of strings to extract multiple
       spectral orders.

    logdamp : int, float, list, tuple
       This gets passed to the optimizer, so its properties are described
       in those routines

    mskorders : many types, optional
       The orders to mask when doing the mult-ended extraction, and this
       can take on different types.  If it is a list, set, or tuple, then
       it is interpreted as multiple orders.  If it is the string 'all',
       then all of the orders in the configuration are used.  If it is
       set as `None`, then no orders are masked.

    algorithm : str, optional
       The optimization algorithm.  Default is 'golden'


    Notes
    -----
    See the paper Ryan, Casertano, & Pirzkal (2018)
    https://ui.adsabs.harvard.edu/abs/2018PASP..130c4501R/abstract
    for more details on the algorithm.  But in brief, the flux in the
    WFSS pixels is modeled as a weighted sum over all sources, all
    wavelengths, and all spectral orders.  Therefore if the weighting
    elements can be determined, then the spectra can be directly
    inferred by inverting the linear system of equations.  In practice,
    these weights are estimated in the `tabulate()` procedure, and
    the equations are cast as a linear operator `Matrix`, which is
    iteratively inverted using sparse-linear algebra techniques.

    """

    DESCRIPTION = "Extracting (multi)"

    def __init__(self, extorders, logdamp, mskorders=None, algorithm='golden',
                 root=None, **kwargs):
        Module.__init__(self, self.extract, **kwargs, multiprocess=False)

        # file root name
        if not isinstance(root, str):
            meta = get_metadata()
            self.root = meta['Name']
        else:
            self.root = root

        self.extorders = as_iterable(extorders)
        self.mskorders = as_iterable(mskorders)
        self.optimizer = optimizer(algorithm, logdamp)
        self.matrix = Matrix(self.extorders, **kwargs)

    def extract(self, data, sources, groups=None, computemodel=True, gzip=True):
        """
        Method to do the mult-ended spectral extraction

        Parameters
        ----------
        data : `su.core.wfss.WFSSCollection`
            The collection of WFSS images

        sources : `su.core.sources.SourceCollection`
            The collection of sources

        groups : `GroupCollection` or None
            The collection of groups.  If set as None, then no grouping
            will be performed.  Default is None

        computemodel : bool
            Flag to compute the model spectral image.  Default is True

        gzip : bool
            Flag to gzip model images.  Default is True

        Notes
        -----
        This is likely not to be directly called.

        """

        # outputting structure
        hdul1 = fits.HDUList()     # simple extractions
        hdul3 = fits.HDUList()     # compound extractions
        hdulL = fits.HDUList()     # L-curve data

        # make a primary header
        phdu = fits.PrimaryHDU()
        headers.add_preamble(phdu.header,
                             filetype=(' ', 'contents of this file'),
                             exttype=('multi', 'method for extraction'))

        headers.add_software_log(phdu.header)              # about software
        data.update_header(phdu.header)                # about WFSS images
        sources.update_header(phdu.header)                 # about the sources
        self.optimizer.update_header(phdu.header)
        Config().update_header(phdu.header)                # global config info
        if groups is not None:
            groups.update_header(phdu.header)

        # add the primary to the file
        hdul1.append(phdu)
        hdul3.append(phdu)

        # get the package meta data
        meta = get_metadata()
        package = meta.get('Name', '')

        # put loops over groups here
        if groups:
            LOGGER.info('Using a group catalog')
        else:
            LOGGER.info('No grouping')
            groups = GroupCollection()

        # container to hold the model values
        models = {}

        # open a PDF to write Grouping images
        pdffile = f'{self.root}_{SUFFIXES["L-curve"]}.pdf'
        LOGGER.info(f'Writing grouped L-curve figure: {pdffile}')
        with PdfPages(pdffile) as pdf:
            # add some info to the PDF
            d = pdf.infodict()
            d['Title'] = 'L-Curve Results'
            d['Author'] = getpass.getuser()
            d['Subject'] = f'L-Curve results for grouped data from {package}'
            d['Keywords'] = f'{package} WFSS L-curve groups'
            d['Producer'] = package

            # keep a list of segids that have not had their spectrum measured.
            # so to start, no objects have been measured, so initialize to
            # the full list of segids
            segids = list(sources.keys())

            # process each group
            for grpid, grpsources in enumerate(groups.groups(sources)):

                # build a matrix for these sources
                self.matrix.build_matrix(data, grpsources, group=grpid)

                # only continue if the matrix is valid
                if self.matrix:

                    # solve the system
                    res = self.optimizer(self.matrix)
                    unc = self.matrix.compute_uncertainty()

                    # initialize a counter
                    s = slice(0, 0)

                    # do each source
                    for segid, source in grpsources.items():
                        # get the extration properties for this source
                        if hasattr(source, 'extpars'):
                            wave = source.extpars.wavelengths()
                        else:
                            wave = self.matrix.defpars.wavelengths()
                        nwave = len(wave)

                        segids.remove(segid)

                        # do each spectral region
                        for regid, region in enumerate(source):
                            s = slice(s.stop, s.stop + nwave)

                            # get the fluxes from the objects
                            flam = res.x[s]
                            func = unc[s]

                            # if something has unc == 0, then it's bad
                            uncons = self.matrix.unconstrained[s]
                            flam[uncons] = np.nan
                            func[uncons] = np.nan

                            # update the outputting structures
                            sources[segid].grpid = grpid
                            sources[segid][regid].sed.reset(wave, flam,
                                                            func=func)

                    # compute the model values
                    if computemodel:
                        model = self.matrix.compute_model(res.x)

                        for key, mod in model.items():
                            if key not in models:
                                models[key] = []
                            models[key].extend(mod)

                # update results for the L-curve data
                kwargs = {'nobj': (len(grpsources), 'Number of sources')}
                hdulL.append(self.matrix.lcurve.as_HDU(grpid=grpid, **kwargs))
                self.matrix.lcurve.pdfplot(pdf, grpid=grpid)

        if computemodel:
            # issue a logging message
            LOGGER.info("Making model images")

            # now let's compute the residuals from the linear model
            for datum in tqdm(data, desc='Rendering models', total=len(data),
                              dynamic_ncols=True):

                # get the output image
                modfile = f'{datum.dataset}_mod.fits'

                # for output file
                hdul = fits.HDUList()
                hdul.append(fits.PrimaryHDU(header=datum.primaryheader()))

                # process each extension
                for detname, detdata in datum.items():
                    # get some info
                    hdr = detdata.headfits('science')
                    unc = detdata.readfits('uncertainty')
                    mod = np.zeros_like(unc, dtype=float)

                    # sum over the model, but recall we need the uncertainty
                    # to undo the inverse-variance weighting ineherent in
                    # the matrix inversion algorithm (see Ryan+ 2018, PASP)
                    key = (datum.dataset, detname)
                    for y, x, m in models[key]:
                        mod[y, x] += m * unc[y, x]

                    # udpate the header info
                    hdr['EXTNAME'] = 'MOD'

                    # save to the output
                    hdul.append(fits.ImageHDU(data=mod, header=hdr))

                # write the file to disk
                hdul.writeto(modfile, overwrite=True)

                # gzip the file if asked
                if gzip:
                    compression.compress(modfile)

        # remove the sources for which we didn't measure a spectrum
        for segid in segids:
            if hasattr(source, 'extpars'):
                wave = source.extpars.wavelengths()
            else:
                wave = self.matrix.defpars.wavelengths()

            nans = np.full_like(wave, np.nan)

            # update the outputting structures
            sources[segid].grpid = -1
            sources[segid][regid].sed.reset(wave, nans, func=nans)

        # loop over sources for outputting
        for source in sources.values():

            hdu = source.as_HDU()
            if source.is_compound:
                hdul3.append(hdu)
            else:
                hdul1.append(hdu)

        # write out the files
        if len(hdul1) > 1:
            hdul1[0].header['FILETYPE'] = '1d spectra'

            x1dfile = f'{self.root}_{SUFFIXES["1d spectra"]}.fits'
            LOGGER.info(f'Writing 1d extractions: {x1dfile}')
            hdul1.writeto(x1dfile, overwrite=True)
        else:
            LOGGER.warning('No 1d spectra written')

        if len(hdul3) > 1:
            hdul3[0].header['FILETYPE'] = '3d spectra'

            x3dfile = f'{self.root}_{SUFFIXES["3d spectra"]}.fits'
            LOGGER.info(f'Writing 3d extractions: {x3dfile}')
            hdul1.writeto(x3dfile, overwrite=True)
        else:
            LOGGER.warning('No 3d spectra written')

        lcvfile = f'{self.root}_{SUFFIXES["L-curve"]}.fits'
        LOGGER.info(f'Writing L-curve tabular data {lcvfile}')
        hdulL.writeto(lcvfile, overwrite=True)
