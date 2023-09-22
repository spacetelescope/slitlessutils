from astropy.io import fits

from ....logger import LOGGER
from . import utils


def _downgrade_wcs_one(imgfile, key, mode, newfile, inplace):
    if not isinstance(imgfile, str):
        LOGGER.warning('imgfile must be a string')
        return None

    LOGGER.info(f'downgrading WCS in {imgfile} to WCSNAME{key}')
    with fits.open(imgfile, mode=mode) as hdul:

        # do some error checking... On Second thought, we probably don't
        # need/want this.
        # obstype = hdul[0].header['OBSTYPE']
        # if obstype != 'IMAGING':
        #    LOGGER.warning(f'direct image obstype ({obstype}) is not imaging.')
        #    return

        for ext, hdu in enumerate(hdul):

            wcsname = hdu.header.get('WCSNAME')

            # do we check for Gaia?  Prolly not.
            # if wcsname and (wcsname.lower().find('gaia') != -1):

            # get the current WCS
            crval = utils.get_crval(hdu.header, '')
            cd = utils.get_cd(hdu.header, '')
            # crvals = [hdu.header['CRVAL1'], hdu.header['CRVAL2']]
            # cd = [[hdu.header['CD1_1'], hdu.header['CD1_2']],
            #       [hdu.header['CD2_1'], hdu.header['CD2_2']]]

            # get the old WCS
            wcsname2 = hdu.header[f'WCSNAME{key}']
            crval2 = utils.get_crval(hdu.header, key)
            cd2 = utils.get_cd(hdu.header, key)
            # crvals2 = [hdu.header[f'CRVAL1{key}'], hdu.header[f'CRVAL2{key}']]
            # cd2 = [[hdu.header[f'CD1_1{key}'], hdu.header[f'CD1_2{key}']],
            #        [hdu.header[f'CD2_1{key}'], hdu.header[f'CD2_2{key}']]]

            # swap the WCSs
            hdul[ext].header['WCSNAME'] = wcsname2
            utils.set_crval(hdul[ext].header, crval2)
            utils.set_cd(hdul[ext].header, cd2)

            # hdul[ext].header['CRVAL1'] = crvals2[0]
            # hdul[ext].header['CRVAL2'] = crvals2[1]
            # hdul[ext].header['CD1_1'] = cd2[0][0]
            # hdul[ext].header['CD1_2'] = cd2[0][1]
            # hdul[ext].header['CD2_1'] = cd2[1][0]
            # hdul[ext].header['CD2_2'] = cd2[1][1]

            hdul[ext].header['WCSNAME'] = wcsname
            utils.set_crval(hdul[ext].header, crval)
            utils.set_cd(hdul[ext].header, cd)

            # hdul[ext].header[f'WCSNAME{key}'] = wcsname
            # hdul[ext].header[f'CRVAL1{key}'] = crvals[0]
            # hdul[ext].header[f'CRVAL2{key}'] = crvals[1]
            # hdul[ext].header[f'CD1_1{key}'] = cd[0][0]
            # hdul[ext].header[f'CD1_2{key}'] = cd[0][1]
            # hdul[ext].header[f'CD2_1{key}'] = cd[1][0]
            # hdul[ext].header[f'CD2_2{key}'] = cd[1][1]

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
        return _downgrade_wcs_one(inputs, key, mode, newfile, inplace)
    elif isinstance(inputs, (list, tuple)):
        return [_downgrade_wcs_one(i, key, mode, newfile, inplace) for i in inputs]
    else:
        LOGGER.warning(f"Invalid datatype ({type(inputs)}) for input.")
        return
