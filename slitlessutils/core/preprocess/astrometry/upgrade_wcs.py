from astropy.io import fits
from astropy.wcs import WCS
import numpy as np

from ....logger import LOGGER
from ...utilities import headers
from ...wfss import WFSS

def upgrade_wcs(imgfile,wfssfile,key='A',newfile=None,inplace=False):
    """
    Method to upgrade the WCS in a WFSS to match that of a direct image 
    that has been tweaked to match Gaia.

    This is necessary for HST *ONLY* as at some point the WCS in the
    direct image was updated with respect to Gaia, but the WFSS images
    were not.  Therefore, to make them match, this method computes the
    multiplicative matrix one needs to transform from the pre-Gaia
    (labeled with subscript 'Gaia') CD matrix to the post-Gaia (labeled
    with the subscript 'IDC' matrix:

    .. math::
       A = CD_{IDC} CD_{Gaia}^{-1}

    This is then matrix-multiplied into the CD matrix for the WFSS image(s).
    In general, this matix A is not simply a rotation matrix as the
    astrometric tweaks applied to the IDC-based astrometry may have included
    skew terms.

    In a similar fasion, the difference in the CRVALs from the direct
    image(s) are used to update the CRVALs in the WFSS image(s):

    .. math::
      \\Delta CRVAL = CRVAL_{Gaia}- CRVAL_{IDC}

    In general, these tweaks are typically <~0.1 arcsec.


    Parameters
    ----------
    imgfile : str
        Full path to the direct image for computing WCS tweaks

    wfssfile : str
        Full path to the WFSS image for applying the WCS tweaks

    key : str, optional
        The WCS key for the pre-Gaia (ie. IDC-based) astrometry. Default
        is 'A'.

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

    LOGGER.info(f'Upgrading WCS in WFSS ({wfssfile}) from direct ({imgfile})')
    if inplace:
        mode = 'update'
    else:
        mode = 'readonly'

    wfss=WFSS.observed(wfssfile)
    with fits.open(imgfile,mode='readonly') as dhdul,\
         fits.open(wfssfile,mode=mode) as whdul:

        # some quick error checking
        if dhdul[0].header['OBSTYPE']!='IMAGING':
            LOGGER.warning('direct image obstype is not imaging.')
            return
        if whdul[0].header['OBSTYPE']!='SPECTROSCOPIC':
            LOGGER.warning('WFSS image obstype is not spectroscopic.')
            return

        # process all the data extensions        
        for (detname,extname),extdata in wfss.extensions():
            exten=extdata.extension
            hdr=dhdul[exten].header
            wcsname=hdr.get(f'WCSNAME{key}')
            
            if wcsname and (wcsname.lower().find('gaia')!=-1):
                # get the best WCS and old WCS
                bwcs = WCS(hdr, dhdul, relax=True)
                owcs = WCS(hdr, dhdul, relax=True, key=key)

                # compute matrix product:
                # CD_best = A . CD_old
                # A = CD_best . inv(CD_old)
                cinv = np.linalg.inv(owcs.wcs.cd)
                A = np.dot(bwcs.wcs.cd, cinv)

                # find shift
                dcrval = (bwcs.wcs.crval-owcs.wcs.crval)

                # get the old WFSS WCS and update it
                wwcs = WCS(whdul[exten].header, whdul, relax=True, key=key)

                # compute corrected WCS as
                # CD_cor = A . CD_WFSS
                # CRVAL_cor = CRVAL_WFSS + (CRVAL_best - CRVAL_old)
                wwcs.wcs.cd = np.dot(A, wwcs.wcs.cd)
                wwcs.wcs.crval = wwcs.wcs.crval+dcrval

                # update the header
                whdul[exten].header['CRVAL1'] = wwcs.wcs.crval[0]
                whdul[exten].header['CRVAL2'] = wwcs.wcs.crval[1]
                whdul[exten].header['CD1_1'] = wwcs.wcs.cd[0, 0]
                whdul[exten].header['CD1_2'] = wwcs.wcs.cd[0, 1]
                whdul[exten].header['CD2_1'] = wwcs.wcs.cd[1, 0]
                whdul[exten].header['CD2_2'] = wwcs.wcs.cd[1, 1]

                # update with some comments
                comment = 'tweak w.r.t. Gaia'
                whdul[exten].header.set('WCSTWEAK', value=True,
                                        comment='Was WCS tweaked w.r.t. Gaia')
                whdul[exten].header.set('WCSTYPE', value='upgrade',
                                        comment = 'WCS upgraded to match Gaia')
                whdul[exten].header.set('DCRVAL1', value=dcrval[0],
                                        comment=comment)
                whdul[exten].header.set('DCRVAL2', value=dcrval[1],
                                        comment=comment)
                whdul[exten].header.set('A1_1', value=A[0, 0], comment=comment)
                whdul[exten].header.set('A1_2', value=A[0, 1], comment=comment)
                whdul[exten].header.set('A2_1', value=A[1, 0], comment=comment)
                whdul[exten].header.set('A2_2', value=A[1, 1], comment=comment)
                headers.add_stanza(whdul[exten].header, 'Tweaked WCS',
                                   before='WCSTWEAK')

        # figure out how to write the file and get a variable to return
        if inplace:
            outfile = wfssfile
        else:
            if newfile is None:
                # get a new filename
                base = os.path.splitext(os.path.basename(wfssfile))[0]
                newfile = f'{base}_twcs.fits'
            whdul.writeto(newfile, overwrite=True)
            outfile = newfile
    return outfile
