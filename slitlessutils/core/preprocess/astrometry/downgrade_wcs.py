from astropy.io import fits

from ....logger import LOGGER
from . import utils


def _downgrade_wcs(imgfile, key, mode, newfile, inplace):
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

    LOGGER.info(f'downgrading WCS in {imgfile} to WCSNAME{key}')
    with fits.open(imgfile, mode=mode) as hdul:

        for ext, hdu in enumerate(hdul):

            wcsname = hdu.header.get('WCSNAME')

            # get the current WCS
            crval = utils.get_crval(hdu.header, '')
            cd = utils.get_cd(hdu.header, '')

            # get the old WCS
            wcsname2 = hdu.header[f'WCSNAME{key}']
            crval2 = utils.get_crval(hdu.header, key)
            cd2 = utils.get_cd(hdu.header, key)

            # swap the WCSs
            hdul[ext].header['WCSNAME'] = wcsname2
            utils.set_crval(hdul[ext].header, crval2)
            utils.set_cd(hdul[ext].header, cd2)

            hdul[ext].header['WCSNAME'] = wcsname
            utils.set_crval(hdul[ext].header, crval)
            utils.set_cd(hdul[ext].header, cd)

            # update with some comments
            hdul[ext].header.set('WCSTWEAK', value=True,
                                 comment='WCS tweaked by slitlessutils')
            hdul[ext].header.set('WCSTYPE', value='downgrade')

        hdul[0].header.add_history('Downgraded astrometry')

        # get a new file
        outfile = utils.new_filename(imgfile, newfile=newfile, inplace=inplace,
                                     suffix=f'WCSNAME{key}')
        if not inplace:
            hdul.writeto(outfile, overwrite=True)

    return outfile


def downgrade_wcs(inputs, key='A', newfile=None, inplace=False):
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

    Returns
    -------
    outfile : str
        The name of the tweaked file.

    """

    mode = 'update' if inplace else 'readonly'

    if isinstance(inputs, str):
        return _downgrade_wcs(inputs, key, mode, newfile, inplace)
    elif isinstance(inputs, (list, tuple)):
        return [_downgrade_wcs(i, key, mode, newfile, inplace) for i in inputs]
    else:
        LOGGER.warning(f"Invalid datatype ({type(inputs)}) for input.")
        return
