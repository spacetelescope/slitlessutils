import os

from astropy.io import fits
from astropy.stats import SigmaClip
import numpy as np

from .....config import SUFFIXES, Config
from .....info import __code__
from .....logger import LOGGER
from ....tables import PDTFile
from ....utilities import headers
from ...module import Module
from .contamination import Contamination

from .spectraltable import SpectralTable
from .boxcar import boxcar


class Single(Module):
    """
    Module to implement single-ended extraction, which is capable of
    the Horne (1986) optimal approach.

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
        A path where output products will be written.  If set to None or
        invalid, then CWD will be used.  If set to valid str, then path
        will be created (if possible), then used.  If still not valid, then
        will use CWD.  Default is None.

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
            self.clipper = SigmaClip(sigma=nsigmaclip, maxiters=maxiters)
        else:
            self.clipper = None

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
        spectra : dict
            A dictionary of `SpectralTable` objects, one for each source

        data : `su.core.wfss.WFSSCollection`
            The collection of WFSS images

        sources : `su.core.sources.SourceColection`
            The collection of sources

        Notes
        -----
        This is likely not to be directly called.
        """

        LOGGER.info("Combining the 1d spectra")

        # stuff for the output data
        prekwargs = {'filetype': (self.FILETYPE, 'contents of this file'),
                     'exttype': ('single', 'method for extraction'),
                     'savecont': (self.savecont, 'Was contamination saved?')}

        # grab the default extraction parameters
        defpars = data.get_parameters()

        # get the flux funits
        config = Config()

        # make an output structure
        hdul = fits.HDUList()

        # make a primary header and fill it with useful info
        phdu = fits.PrimaryHDU()
        headers.add_preamble(phdu.header, **prekwargs)  # basic preamble
        headers.add_software_log(phdu.header)           # software props
        data.update_header(phdu.header)                 # WFSS image props
        sources.update_header(phdu.header)              # source props
        config.update_header(phdu.header)               # global config info
        self.contamination.update_header(phdu.header)   # contamination prop

        # if there is a sigma clipper
        phdu.header.set('SIGCLIP', value=bool(self.clipper),
                        comment='Sigma-clip weighting?')
        if phdu.header['SIGCLIP']:
            phdu.header.set('NSIGMA', value=self.clipper.sigma,
                            comment='number of sigma for clipping')
            phdu.header.set('MAXITER', value=self.clipper.maxiters,
                            comment='maximum number of iterations')
        headers.add_stanza(phdu.header, 'Spectral Combination',
                           before='SIGCLIP')
        super().update_header(phdu.header)            # about the CPU settings

        # add to the HDUList
        hdul.append(phdu)

        # now process each object
        for segid, source in sources.items():

            # set the spectral extraction parameters
            pars = defpars.update_pars(source[0])

            # get all the spectra for this object
            spectra = SpectralTable(segid)
            for result in results:
                if segid in result:
                    spectra += result[segid]

            # combine it
            spectrum = spectra.combine(pars, clipper=self.clipper)

            # update the source
            source[0].sed.reset(spectrum.get('wave'),
                                spectrum.get('flam'),
                                func=spectrum.get('func'),
                                cont=spectrum.get('cont'),
                                npix=spectrum.get('npix'))

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
                out = tuple(np.empty(0, dtype=a.dtype) for a in args)
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
        spectra : `SpectralTable`
            The measurements for all the spectra for this source and is passed to `self._combine`

        Notes
        -----
        This is likely not to be directly called.

        """

        # get the defaults

        # sort out the extraction mode
        extmode = kwargs.get('extmode', 'boxcar').lower()

        # get the profile for horne?
        # profile = kwargs.get('profile', 'uniform')

        # padding for the contamination models
        padx = kwargs.get('padx', (5, 5))
        pady = kwargs.get('pady', (5, 5))
        width = kwargs.get('width', 5)

        # create a data structure to hold the results
        spectra = {segid: SpectralTable(segid) for segid in sources.keys()}

        # contamination image
        if self.savecont:
            chdul = fits.HDUList()
            phdr = fits.PrimaryHDU()

            LOGGER.debug("add stuff to contam images' primary headers")
            chdul.append(phdr)

        # open the PDT for reading
        with PDTFile(data, path=self.path, mode='r') as h5:

            for detname, detdata in data.items():
                # load the tables and configuration
                h5.load_detector(detname)
                cfg = detdata.config

                # load the data
                sci, unc, dqa = detdata.readimages()

                # get the flat field
                flatfield = detdata.config.load_flatfield(**kwargs)

                # initialize the contamination modeling
                self.contamination.make_model(sources, h5)

                # process each order
                orders = self.extorders if self.extorders else detdata.orders
                for ordname in orders:
                    h5.load_order(ordname)
                    order = cfg[ordname]

                    # process each source
                    for segid, source in sources.items():

                        odt = h5.load_odt(source)
                        if odt:
                            # get the extraction parameters
                            # pars = defpars.update_pars(source[0])

                            # do a contmaination if requestied
                            if self.contamination:
                                # get bounding box
                                x0, x1, y0, y1 = odt.bounding_box(dx=padx,
                                                                  dy=pady)

                                # get the contam data
                                chdu = self.contamination(segid, ordname,
                                                          h5, sources,
                                                          detdata, flatfield,
                                                          bbx=(x0, x1), bby=(y0, y1))

                                if self.savecont:
                                    chdul.append(chdu)

                                # get the offsets in the contamination image
                                # xoff = -chdu.header.get('LTV1', 0)
                                # yoff = -chdu.header.get('LTV2', 0)

                            if extmode == 'boxcar':
                                spectrum = boxcar(source, detdata, sci, unc, dqa, flatfield, order,
                                                  odt, width=width)

                                spectra[segid] += spectrum

                            elif extmode == 'horne':
                                raise NotImplementedError

        if self.savecont:
            chdul.writeto(f'{data.dataset}_cont.fits', overwrite=True)

        return spectra
