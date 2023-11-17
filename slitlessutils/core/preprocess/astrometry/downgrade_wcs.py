from astropy.io import fits

from ....logger import LOGGER
from . import utils


def _downgrade_wcs(imgfile, key, mode, newfile, inplace, force):
    """
    Helper function to actually execute the downgrading

    Parameters
    ----------
    imgfile : str
        Path to a file to downgrade

    key : str
        WCS key to downgrade from.  Will replace the current solutions

    mode : str
        The file mode see `astropy.io.fits.open()`

    newfile : str
        The output filename

    inplace : bool
        Flag that the update should happen in place

    Returns
    -------
    outfile : str
        The name of the file that contains the new WCS

    Notes
    -----
    It is not expected for one to directly call this function, though nothing prohibits it.

    """

    if not isinstance(imgfile, str):
        LOGGER.warning('imgfile must be a string')
        return None

    with fits.open(imgfile, mode=mode) as hdul:
        # check if already tweaked
        wcscorr = hdul[0].header.get('WCSCORR', False)
        if wcscorr and not force:
            msg = f'Astrometry was already corrected: {imgfile}'
            LOGGER.warning(msg)
            return imgfile
        LOGGER.info(f'downgrading WCS in {imgfile} to WCSNAME{key}')

        for ext, hdu in enumerate(hdul):

            extname = hdu.header.get('EXTNAME')
            if extname == 'SCI':

                # get the current WCS
                wcsname = hdu.header.get('WCSNAME')
                crval = utils.get_crval(hdu.header, '')
                cd = utils.get_cd(hdu.header, '')

                # get the old WCS
                wcsname2 = hdu.header.get(f'WCSNAME{key}')
                crval2 = utils.get_crval(hdu.header, key)
                cd2 = utils.get_cd(hdu.header, key)

                # swap the WCSs
                hdul[ext].header['WCSNAME'] = wcsname2
                utils.set_crval(hdul[ext].header, crval2, '')
                utils.set_cd(hdul[ext].header, cd2, '')

                hdul[ext].header[f'WCSNAME{key}'] = wcsname
                utils.set_crval(hdul[ext].header, crval, key)
                utils.set_cd(hdul[ext].header, cd, key)

                # update with some comments
                utils.update_wcshistory(hdul[ext].header, 'downgrade')

        utils.update_wcshistory(hdul[0].header, 'downgrade')

        # get a new file
        outfile = utils.new_filename(imgfile, newfile=newfile, inplace=inplace,
                                     suffix=f'WCSNAME{key}')
        if not inplace:
            hdul.writeto(outfile, overwrite=True)

    return outfile


def downgrade_wcs(inputs, key='A', newfile=None, inplace=False, force=False):
    """
    Method to take WCS in a direct image that has been tied to Gaia, and
    downgrade it to a previous WCS.

    Parameters
    ----------
    inputs : str, list, or tuple
        If this is a str, then it is the name of the file to be downgrade.
        If it is a list or tuple, then it is many such filenames.

    key : str, optional
        the WCS key for the old WCS to downgrade to. default is 'A'

    inplace : bool, optional
        A flag that the tweaks should be applied to the WFSS image in place,
        ie. without creating a new file.  Default is False

    newfile : str, optional
        See the file `utils.py` for more details.  Default is None

    force : bool, optional
        Force the updating of the astrometry, even if done already.
        Default is False.

    Returns
    -------
    outfile : str
        The name of the tweaked file.

    """

    mode = 'update' if inplace else 'readonly'

    if isinstance(inputs, str):
        return _downgrade_wcs(inputs, key, mode, newfile, inplace, force)
    elif isinstance(inputs, (list, tuple)):
        return [_downgrade_wcs(i, key, mode, newfile, inplace, force) for i in inputs]
    else:
        LOGGER.warning(f"Invalid datatype ({type(inputs)}) for input.")
        return
