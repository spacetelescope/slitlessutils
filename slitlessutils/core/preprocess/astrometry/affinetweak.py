import numpy as np
from astropy.io import fits

from ....logger import LOGGER
from . import utils


class AffineTweak(dict):
    """
    Class to contain an affine perturbation to the Astrometry of a fits image
    between two different WCS solutions present in a file.

    This assumes the "classic" astrometric formulation of a CD (or PC) matrix
    and a CRVAL vector.  The affine transformation is a simple vector difference
    of the CRVAL keywords and a matrix product between the CD matrices.  This
    permits for skew, rotation, and translation between the two WCSs.  Now, the
    method computes the multiplicative matrix of the CD-matrix and vector
    difference between the final WCS and the initial WCS:

    .. math::
       A = CD_{0} CD_{1}^{-1}

    This is then matrix-multiplied into the CD matrix for a different image to
    update it.  In general, this matrix A is not simply a rotation matrix as
    the astrometric tweaks may include skew or pixel-size changes.  In a similar
    fashion, the difference in the CRVALs between the final and initial WCSs:

    .. math::
      \\Delta CRVAL = CRVAL_{1} - CRVAL_{0}

    In general, these tweaks are typically <~0.1 arcsec.


    Parameters
    ----------
    reffile : str
        The file name of the image to use as a reference.

    key0 : str, optional
        The WCS key for the initial WCS.  This is usually the "older" WCS.
        Note, this is not case-sensitive (only capitals are used).
        Default is "A"

    key1 : str, optional
        The WCS key for the final WCS.  This is usually the "newer" WCS.
        Note, this is not case-sensitive (only capitals are used).
        Default is '' (the empty string).

    Notes
    -----
    This class inherits from ``dict`` to store individual tweaks for each
    detector, in the case of instruments that have multiple detectors (e.g.
    ACS WFC).  The keyword/value pairs for the dict will be the extension
    name/version and a dict of data for the perturbations.  Users are
    highly discouraged from manually manipulating those contents.

    Examples
    --------
    The usual work flow would be to pass a direct image that has been
    tweaked to some optimal solution (e.g. referenced to Gaia) to derive
    the perturbations.  Then apply those perturbations to a grism image
    that has not been similarly updated.

    >>> tweak = AffineTweak('directimage_flt.fits') # doctest: +SKIP
    >>> newfile = tweak('grismimage_flt.fits') # doctest: +SKIP

    """

    def __init__(self, reffile, key0='A', key1=''):

        self.reffile = reffile
        self.key0 = key0.upper()
        self.key1 = key1.upper()

        with fits.open(self.reffile, mode='readonly') as hdul:
            for hdu in hdul:
                # get the extension name
                extname = hdu.header.get('EXTNAME')
                extver = hdu.header.get('EXTVER')
                ext = (extname, extver)

                if extname == 'SCI':
                    # get the old WCS
                    cd0 = utils.get_cd(hdu.header, self.key0)
                    crval0 = utils.get_crval(hdu.header, self.key0)

                    # get the new WCS
                    cd1 = utils.get_cd(hdu.header, self.key1)
                    crval1 = utils.get_crval(hdu.header, self.key1)

                    # create and fill a dict of the data
                    self[ext] = {}
                    self[ext]['A'] = np.dot(cd1, np.linalg.inv(cd0))
                    self[ext]['d'] = crval1 - crval0
                    self[ext]['wcsname0'] = hdu.header[f'WCSNAME{self.key0}']
                    self[ext]['wcsname1'] = hdu.header[f'WCSNAME{self.key1}']

    def __call__(self, datfile, inplace=False, newfile=None, force=False):
        """
        Main method to apply the tweaks to a given file.

        Parameters
        ----------
        datfile : str
            The name of the file to be tweaked.  This should have the same
            number and names of the extensions as the 'reffile' used to
            instantiate the object.

        inplace : bool, optional
            A flag that this file should be updated in place.  Note, if
            the WCS in a file is overwritten (ie. `inplace=True`), then
            the original WCS is unrecoverable.  Default is False.

        newfile : str or None, optional
            See the file `utils.py` for more details.  Default is None.

        force : bool, optional
            Force the updating of the astrometry, even if done already.
            Default is False.

        Returns
        -------
        filename : str
            The name of the output file that was updated.

        """

        mode = 'update' if inplace else 'readonly'
        with fits.open(datfile, mode=mode) as hdul:
            # check if already tweaked
            wcscorr = hdul[0].header.get('WCSCORR', False)
            if wcscorr and not force:
                msg = f'Astrometry was already corrected: {datfile}'
                LOGGER.warning(msg)
                return

            LOGGER.info(f'upgrading WCS in {datfile} from {self.reffile}')
            for ext, data in self.items():

                wcsname0 = hdul[ext].header[f'WCSNAME{self.key0}']
                if wcsname0 != data['wcsname0']:
                    LOGGER.warning('mismatches in the WCSNAMEs')

                cd = utils.get_cd(hdul[ext].header, self.key0)
                crval = utils.get_crval(hdul[ext].header, self.key0)

                cd = np.dot(data['A'], cd)
                crval = crval + data['d']

                hdul[ext].header.set('WCSNAME', data['wcsname1'])
                utils.set_cd(hdul[ext].header, cd, self.key1)
                utils.set_crval(hdul[ext].header, crval, self.key1)

                # add the tweak to the header
                utils.update_wcshistory(hdul[ext].header, 'affine upgrade')
                hdul[ext].header.set('DCRVAL1', value=data['d'][0],
                                     comment='Change in CRVAL1 (deg)')
                hdul[ext].header.set('DCRVAL2', value=data['d'][1],
                                     comment='Change in CRVAL2 (deg)')
                hdul[ext].header.set('A1_1', value=data['A'][0, 0],
                                     comment='Affine transform')
                hdul[ext].header.set('A1_2', value=data['A'][0, 1],
                                     comment='Affine transform')
                hdul[ext].header.set('A2_1', value=data['A'][1, 0],
                                     comment='Affine transform')
                hdul[ext].header.set('A2_2', value=data['A'][1, 1],
                                     comment='Affine transform')

            utils.update_wcshistory(hdul[0].header, 'affine upgrade')

            outfile = utils.new_filename(datfile, inplace=inplace,
                                         suffix='twcs')
            if not inplace:
                hdul.writeto(outfile, overwrite=True)

            return outfile
