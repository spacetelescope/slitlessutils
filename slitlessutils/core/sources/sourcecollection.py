import numpy as np
from astropy.wcs import WCS, FITSFixedWarning, utils as wcsutils
from astropy.io import fits


from contextlib import nullcontext
import os
import warnings


from .sedfile import SEDFile
from .source import Source
from ..photometry import Throughput
from ..utilities import indices, headers
from ...logger import LOGGER


class SourceCollection(dict):
    """
    Class to hold a bunch of sources.

    inherits from dict, where the key/value are the segid and the source
    object

    """

    # when cropping a source, this is the extra padding.  This should be
    # kept somewhat large so that a local sky background can be reliably
    # computed.  this is in units of pixels
    PAD = 20

    # the default zeropoint
    DEFZERO = 26.

    def __init__(self, segfile, detfile, maglim=np.inf, minpix=0, preprocessor=None,
                 zeropoint=None, throughput=None, sedfile=None,
                 **kwargs):
        """
        Initializer

        Parameters
        ----------
        segfile : str
            Full path to the segmentation image

        detfile : str
            Full path to the direct image

        maglim : float or None, optional
            The AB magnitude limit to eliminate faint sources.
            Default is `np.inf`.

        minpix : int, optional
            The minimum number of pixels for a segmentation region to be
            considered valid.  Default is 0

        preprocessor : callable, optional
            A preprocessing function.  see `su.core.sources.operators.py`
            for more information.  Default is None

        zeropoint : int, float, optional
            The AB magnitude zeropoint for the detection image.  If `None`
            then the header of the detection image is searched.  The
            possible keywords are:

            ('ABZERO','ZEROPT','MAGZERO','ZERO','ZPT')

            See the method: `get_zeropoint()` for more info.
            If the telescope is JWST, then a zeropoint is computing
            based on header values (since those images are in MJy/sr)
            Default is 26.

        throughput : `su.core.photometry.Throughput` or None, optional
            The throughput object.  If `None`, then it is looked for
            based on header keywords (see note below).  Default is None

        sedfile : str, optional
            The full path to a mult-extension fits file of binary
            tabulated data.  This will enable loading SEDs directly
            into the sources, which is primarily used for simulation
            Default is None

        kwargs : dict, optional
            Used to pass many of inputs around


        Notes
        -----
        1) The detection and segmentation images must have the same
           dimensionality and number of extensions

        2) The detection image should have a few header keywords set:
           a) TELESCOP= the telescope name (shall be 'HST' or 'JWST'),
                        but the capitalization does not matter
           b) INSTRUME= the instrument to use.
           c) FILTER  = the filter to use.  This should be a string
                        such as 'F125W'


           -> these three keywords are used to load a `Throughput` object
              with the classmethod: `Throughput.from_keys()`

           -> the filters are stored in the configuration directory:
              slitlessutils_config/bandpasses/


        """

        self.segfile = segfile
        self.detfile = detfile
        self.maglim = maglim
        self.minpix = minpix
        self.preprocessor = preprocessor
        self.sedfile = sedfile

        fitsargs = {'mode': 'readonly', 'ignore_missing_simple': True}
        with fits.open(self.segfile, **fitsargs) as hdus,\
                fits.open(self.detfile, **fitsargs) as hdud:

            # hide warnings that are an issue for JWST
            warnings.simplefilter("ignore", FITSFixedWarning)

            # sort out the input nature
            nhdus, nhdud = len(hdus), len(hdud)
            if nhdus != nhdud:
                msg = f'Segmentation {self.segfile} and direct ' +\
                    f'{self.detfile} have differing number of ' +\
                    f'extensions: {nhdus},{nhdud}'
                LOGGER.error(msg)
                return

            # get the ZEROPOINT and THROUGHPUT
            self.zeropoint = self.get_zeropoint(hdud, zeropoint=zeropoint)
            self.throughput = self.get_throughput(hdud, throughput=throughput, **kwargs)

            # is it an MEF?
            self.mef = nhdus > 1
            if self.mef:
                self._load_mef(hdus, hdud, **kwargs)
            else:
                self._load_classic(hdus, hdud, exten=0, **kwargs)

        # load SEDs if asked
        if sedfile:
            self.load_sedlib(sedfile)

    def _load_classic(self, hdus, hdui, exten=0, **kwargs):
        """
        Method to load a classic segmentation map

        This is unlikely to be directly called.

        Parameters
        ----------
        hdus : `astropy.io.fits.HDUList`
             A list of HDUs for the segmentation image

        hdui : `astropy.io.fits.HDUList`
             A list of HDUs for the direct image

        exten : int, optional
             Specifier to which extension to use

        kwargs : dict, optional
             The optional keywards

        Notes
        -----
        A classic segmentation map is defined as the segmentation image
        produced by default by tools like Source Extractor.


        """

        LOGGER.info(f'Loading a classic segmentation map: {self.segfile}')

        seg = hdus[exten].data
        img = hdui[exten].data
        hdr = hdui[exten].header

        assert (seg.shape == img.shape), \
            f'Images have different dimensions: {seg.shape} {img.shape}'

        # just some numbers for easy access later
        ny, nx = img.shape

        # what is datatype of segmentation image.
        # int are treated as normal segmaps
        # flt are treated as complex sources
        isfloat = issubclass(seg.dtype.type, np.floating)
        if isfloat:
            seg, reg = np.divmod(seg, 1)
            seg = seg.astype(int)
            reg = np.floor(reg*self.MAXSPECREG).astype(int)

        # find pixels for each object
        ri = indices.reverse(seg, ignore=(0,))
        for segid, (y, x) in ri.items():
            # get a bounding box
            x0 = np.maximum(np.amin(x)-self.PAD, 0)
            x1 = np.minimum(np.amax(x)+self.PAD+1, nx-1)
            y0 = np.maximum(np.amin(y)-self.PAD, 0)
            y1 = np.minimum(np.amax(y)+self.PAD+1, ny-1)

            # excise the image
            subseg = seg[y0:y1, x0:x1]
            subimg = img[y0:y1, x0:x1]

            # update the header
            subhdr = hdr.copy()
            subhdr['NAXIS1'] = subimg.shape[1]
            subhdr['NAXIS2'] = subimg.shape[0]
            subhdr['CRPIX1'] -= x0
            subhdr['CRPIX2'] -= y0
            subhdr['LTV1'] = hdr.get("LTV1", 0.)-x0
            subhdr['LTV2'] = hdr.get("LTV2", 0.)-y0

            # make a region id
            if segid < 0:
                # make checkerboard reg
                subreg = np.zeros_like(subseg, dtype=int)
                gy, gx = np.where(subseg == segid)
                subreg[gy, gx] = np.arange(1, len(gx)+1, dtype=int)
            elif isfloat:
                # use existing reg
                subreg = reg[y0:y1, x0:x1]
            else:
                subreg = (subseg == segid).astype(int)
            self.addSource(segid, subimg, subseg, subreg, subhdr)

    def _load_mef(self, hdus, hdui, **kwargs):
        """
        Method to load a multi-extension fits (MEF) file for the
        direct and segmentation images

        This is unlikely to be directly called.

        Parameters
        ----------
        hdus : `astropy.io.fits.HDUList`
             A list of HDUs for the segmentation image

        hdui : `astropy.io.fits.HDUList`
             A list of HDUs for the direct image

        kwargs : dict, optional
             The optional keywards

        Notes
        -----
        A classic segmentation map is defined as the segmentation image
        produced by default by tools like Source Extractor.


        """

        LOGGER.info(f'Loading a MEF segmentation map: {self.segfile}')
        raise NotImplementedError()

    def addSource(self, segid, img, seg, reg, hdr):
        """
        Method to introduce a new source to the collection

        Parameters
        ----------
        segid : int
            The object ID number

        img : `np.ndarray`
            The direct image

        seg : `np.ndarray`
            The segmentation map

        hdr : `astropy.io.fits.Header`
            The fits header for the direct image

        """
        if callable(self.preprocessor):
            img, seg, hdr = self.preprocessor(img, seg, hdr)

        source = Source(segid, img, seg, reg, hdr, zeropoint=self.zeropoint)
        if source and source.mag < self.maglim and source.npixels > self.minpix:

            # put in flat level to SED
            wave = self.throughput.photplam
            for region in source:
                fnu = source.fnu*np.sum(region.w)
                flam = fnu*(2.99e10/wave)*(1e8/wave)  # /Config().fluxscale
                # region.sed.set_sed([wave],[flam])
                region.sed.append(wave, flam)

            # self[hdr['SEGID']]=source
            self[source.segid] = source

    def set_spectral_parameters(self, **kwargs):
        """
        Method to set the spectral parameters for all of the sources that
        are part of this collection

        Parameters
        ----------
        kwargs : dict, optional
            Dictionary of keywords, can be 'wave0','wave1','dwave','scale'

        """
        for segid, source in self.items():
            source.set_spectral_parameters(**kwargs)

    def get_spectral_parameters(self):
        raise NotImplementedError

    def update_header(self, hdr):
        """
        Method to update an existing fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            The fits header to be updated
        """

        maglim = 'INF' if np.isinf(self.maglim) else self.maglim

        sedfile = self.sedfile if self.sedfile else ' '*8

        hdr.set('SEGFILE', value=self.segfile, comment='segmentation image')
        hdr.set('DETFILE', value=self.detfile,
                comment='direct image for profile weights')
        hdr.set('MEF', value=self.mef, comment='is it a MEF file?')

        if hasattr(self, 'sedstype'):
            if self.sedstype == 'sedfile':
                hdr.set('SEDFILE', value=self.sedfile,
                        comment='file *INPUT* containing SEDs')
            elif self.sedstype == 'imglist':
                hdr.set("IMGLIST", value=self.imglist,
                        comment='list of fits files for broadband SEDs')
            elif self.sedstype == 'sedcat':
                hdr.set("SEDCAT", value=self.sedcat,
                        comment='catalog of broadband photometry')

        hdr.set('NSOURCE', value=len(self), comment='number of sources')
        hdr.set('MAGLIM', value=maglim, comment='magnitude limit of sources')
        headers.add_stanza(hdr, 'Source Properties', before='SEGFILE')

        if self.preprocessor is not None:
            self.preprocessor.update_header(hdr)

    def get_zeropoint(self, hdul, zeropoint=None, exten=0):  # ,**kwargs):
        """
        Helper method to get AB zeropoint from the header


        Parameters
        ----------
        hdul : astropy.io.fits.HDUList
           an HDUList opened from `astropy.io.fits.open()`.  This is used
           to grab the header infor from a particular extension (as passed
           by `exten` keyword).

        zeropoint : float or `None`, optional
           This is a default AB zeropoint to use.  This is likely to not be
           directly used.  Default is 26

        exten : int, str, or tuple, optional
           The extension name.  Default is 0

        Returns
        -------
        zeropoint : float
           The AB magnitude zero point

        Notes
        -----
        1) It is unlikely that this is directly called.
        2) The logic is to first return the zeropoint if the optional
           argument `zeropoint` is a numeric type.  If not, then the
           header is inspected for the `TELESCOP` keyword, and different
           things are done based on that value.  If `TELECSCOP=="JWST"`,
           then the zeropoint is computed from:

        .. math::
          Z = -2.5*np.log10((1e6/3631)*(np.pi/180)**2*pixarea)

           If `TELESCOPE=="HST"`, then additional keywords are queried in
           the following order: "ABZERO", "ZEROPT", "MAGZERO", "ZERO", and
           "ZPT".  If none are found then the default of 26.0 is used and
           a warning is issued.
        """

        if zeropoint:
            return zeropoint

        tel = hdul[exten].header.get('TELESCOP', '')

        if tel == 'JWST':
            w = WCS(hdul[exten].header, relax=True)
            pixarea = wcsutils.proj_plane_pixel_area(w)
            zeropoint = -2.5*np.log10((1e6/3631)*(np.pi/180)**2*pixarea)

        elif tel == 'HST':

            for k in ('ABZERO', 'ZEROPT', 'MAGZERO', 'ZERO', 'ZPT'):
                if k in hdul[0].header:
                    zeropoint = hdul[exten].header[k]
                    break
        else:
            pass

        if not zeropoint:
            msg = f"No zeropoint was loaded, using default: {self.DEFZERO}"
            LOGGER.warn(msg)
            zeropoint = self.DEFZERO

        return zeropoint

    def get_throughput(self, hdul, throughput=None, exten=0, **kwargs):
        """
        Helper method to read a throughput curve based the header info

        Parameters
        ----------
        hdul : `astropy.io.fits.HDUList`
           an HDUList opened from `astropy.io.fits.open()`.  This is used
           to grab the header infor from a particular extension (as passed
           by `exten` keyword).

        throughput : None or `slitlessutls.core.photometry.Throughput`, optioanl
           The default throughput to use.


        Returns
        -------
        thru : `slitlessutls.core.photometry.Throughput`
           the throughput

        Notes
        -----
        1) It is unlikely that this is directly called.
        2) If a valid throughput is passed, then use that.  If not, then
           check the header for the keyword `FILTFILE` as the filter file
           on disk as full path.  If that doesn't exist, then check that
           all of "TELESCOP", "INSTURME" and "FILTER" exist, If so, then
           use those to load a throughput.  If still fails, then return
           `None` and issue a warning.
        """

        if throughput:
            return throughput

        if 'FILTFILE' in hdul[exten].header:
            filtfile = kwargs.get('filtfile', hdul[exten].header['FILTFILE'])
            throughput = Throughput.from_file(filtfile)
        elif all(k in hdul[exten].header for k in ('TELESCOP', 'INSTRUME', 'FILTER')):
            tel = kwargs.get('telescope', hdul[exten].header['TELESCOP']).lower()
            ins = kwargs.get('instrument', hdul[exten].header['INSTRUME']).lower()
            fil = kwargs.get('filter', hdul[exten].header['FILTER']).lower()

            thru = Throughput.from_keys(tel, ins, fil)

        else:
            LOGGER.warning("Unable to load any throughput curve")
            thru = None

        return thru

    def load_sedlib(self, sedfile):
        """
        Method to load SEDs from a library file, which is in the fits
        format.

        Parameters
        ----------
        sedfile : `SEDFile`
            The library file
        """

        LOGGER.info(f'loading SEDs from SEDFile: {sedfile}')
        self.sedstype = 'sedfile'
        self.sedfile = sedfile
        with SEDFile(self.sedfile) as sedlib:
            for source in self.values():
                source.load_sedlib(sedlib, throughput=self.throughput)

    def load_sedimages(self, imglist):
        """
        Method to load the SEDs for the sources as a set of direct
        images

        """

        LOGGER.info(f'loading SEDs from imglist: {imglist}')

        # load bunch of space-based images by pixels

        self.sedstype = 'imglist'
        raise NotImplementedError

    def load_sedphotometry(self, sedcat):
        """
        Method to load SEDs into the objects from a comma-separated
        value (CSV)

        Parameters
        ----------
        sedcat : str
            Full path to an ascii file that contains the photometry

        Notes
        -----
        The file should contain a header row that describes how the file
        is interpreted.  The header should be the first line and
        contain the wavelengths (in A), be delimited by a comma, and
        commented by the '#' symbol.

        The data section should start with the segid followed by all the
        fluxes, and again be comma-separated. Therefore the first entry
        in the header row should be blank.  Finally, the fluxes should
        be in *microJanskies*.

        Example
        -------
        <start_of_file>
        #  ,5400,6500,7780,9000
        4,15.,26.,57.,100.
        6,1.,6.,8.,12.
        10,0.,75.,70.,73.

        <end_of_file>


        """

        # load bunch of single band catalog data
        LOGGER.info(f'loading SEDs from sedcat: {sedcat}')
        self.sedstype = 'sedcat'
        self.sedcat = sedcat
        with open(cat, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#'):

                    tokens = line[1:].strip().split(',')
                    waves = np.array([float(t) for t in tokens], dtype=float)
                    nfilt = len(waves)
                    ncol = nfilt+1
                    s = np.argsort(waves)
                    waves = waves[s]
                else:
                    tokens = line.split(',')
                    if len(tokens) == ncol:
                        segid = int(tokens.pop(0))

                        if segid in self:

                            # convert from uJy to flam
                            fnu = np.array(tokens, dtype=float)[s]
                            flam = fnu*1e-29*(2.998e10/waves)*(1e8/waves)

                            # put this in source.load_seds()
                            self[segid].load_sedarray(waves, flam)

    def write_seds(self, **kwargs):
        """
        Method to write the SEDs to individual files.

        Parameters
        ----------
        kwargs : dict, optional
            Parameters passed to `Source().write_seds()`
        """

        for source in self.values():
            source.write_seds(**kwargs)

    def __str__(self):
        return f'Source Collection with {len(self)} sources'
