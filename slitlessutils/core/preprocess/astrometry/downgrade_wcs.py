import os

from astropy.io import fits

from ....logger import LOGGER


def downgrade_wcs(imgfile, key='A', newfile=None, inplace=False):
    """
    Method to take WCS in a direct image that has been tied to Gaia, and
    downgrade it to a previous WCS

    Parameters
    ----------
    imgfile : str
        the name of the direct image file

    key : str, optional
        the WCS key for the old WCS to downgrade to. default is 'A'

    inplace : bool, optional
        A flag that the tweaks should be applied to the WFSS image in place,
        ie. without creating a new file.  Default is False

    newfile : str, optional
        The name of the new file with the tweaked WCS.  If inplace==False,
        then this is ignored.  If inplace==True, then this name will be used.
        If newfile == None, then a filename will be generated as:

        f'{}_twcs.fits'

    Returns
    -------
    outfile : str
        The name of the tweaked file.

    """

    LOGGER.info(f'downgrading WCS in {imgfile} to WCSNAME{key}')

    mode = 'update' if inplace else 'readonly'
    with fits.open(imgfile, mode=mode) as hdul:
        # do some error checking
        obstype = hdul[0].header['OBSTYPE']
        if obstype != 'IMAGING':
            LOGGER.warning(f'direct image obstype ({obstype}) is not imaging.')
            return

        for ext, hdu in enumerate(hdul):

            wcsname = hdu.header.get('WCSNAME')

            if wcsname and (wcsname.lower().find('gaia') != -1):

                # get the current WCS
                crvals = [hdu.header['CRVAL1'], hdu.header['CRVAL2']]
                cd = [[hdu.header['CD1_1'], hdu.header['CD1_2']],
                      [hdu.header['CD2_1'], hdu.header['CD2_2']]]

                # get the old WCS
                wcsname2 = hdu.header[f'WCSNAME{key}']
                crvals2 = [hdu.header[f'CRVAL1{key}'], hdu.header[f'CRVAL2{key}']]
                cd2 = [[hdu.header[f'CD1_1{key}'], hdu.header[f'CD1_2{key}']],
                       [hdu.header[f'CD2_1{key}'], hdu.header[f'CD2_2{key}']]]

                # swap the WCSs
                hdul[ext].header['WCSNAME'] = wcsname2
                hdul[ext].header['CRVAL1'] = crvals2[0]
                hdul[ext].header['CRVAL2'] = crvals2[1]
                hdul[ext].header['CD1_1'] = cd2[0][0]
                hdul[ext].header['CD1_2'] = cd2[0][1]
                hdul[ext].header['CD2_1'] = cd2[1][0]
                hdul[ext].header['CD2_2'] = cd2[1][1]

                hdul[ext].header[f'WCSNAME{key}'] = wcsname
                hdul[ext].header[f'CRVAL1{key}'] = crvals[0]
                hdul[ext].header[f'CRVAL2{key}'] = crvals[1]
                hdul[ext].header[f'CD1_1{key}'] = cd[0][0]
                hdul[ext].header[f'CD1_2{key}'] = cd[0][1]
                hdul[ext].header[f'CD2_1{key}'] = cd[1][0]
                hdul[ext].header[f'CD2_2{key}'] = cd[1][1]

                # update with some comments
                hdul[ext].header.set('WCSTWEAK', value=True,
                                     comment='WCS tweaked by slitlessutils')
                hdul[ext].header.set('WCSTYPE', value='downgrad',
                                     comment='WCS downgraded to pre-Gaia')

        # figure out how to write the file and get a variable to return
        if inplace:
            outfile = imgfile

        else:
            if newfile is None:
                # get a new filename
                base = os.path.splitext(os.path.basename(imgfile))[0]
                newfile = f'{base}_WCSNAME{key}.fits'
            hdul.writeto(newfile, overwrite=True)
            outfile = newfile
    return outfile
