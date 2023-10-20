import datetime
import os

import numpy as np


def new_filename(imgfile, newfile=None, inplace=False, suffix='wcs'):
    """
    Method to get a new filename following some rules

    If ``inplace==True``, then return the input name; else if
    ``newfile`` is passed, then return that; else parse the input
    filename and return an updated name.

    Parameters
    ----------
    imgfile : str
       The filename to manipulate

    newfile : str, optional
       The name of the output file.  If `None`, then a new filename will
       be created.  Default is None.

    inplace : bool, optional
       Flag indicating if file should be updated in place, on disk.
       Default is False.

    suffix : str, optional
       The suffix for the file.  Default is 'wcs'.

    Returns
    -------
    outfile : str
       The name of the output file.
    """

    if inplace:
        outfile = imgfile
    else:
        if newfile is None:
            base = os.path.splitext(os.path.basename(imgfile))[0]
            newfile = f'{base}_{suffix}.fits'
        outfile = newfile

    return outfile


def update_wcshistory(hdr, wcstype):
    """
    Function to update the WCS history in the header

    Parameters
    ----------
    hdr : `astropy.io.fits.Header()`
       The header to update

    wcstype : str
       The name of the WCS type

    """

    now = datetime.datetime.now()
    hdr.set('WCSTWEAK', value=True, comment='WCS tweaked by slitlessutils')
    hdr.set('WCSTYPE', value=wcstype)
    hdr.add_history(f'WCS updated by slitlessutils on {now.isoformat()}."')


def get_cd(hdr, key):
    """
    Function to retrieve a CD matrix from an `astropy.io.fits.Header()`

    Parameters
    ----------
    hdr : `astropy.io.fits.Header()`
       The header to retrieve from

    key : str
       The WCS key for which WCS to pull from

    Returns
    -------
    cd : `np.ndarray()`
       The CD matrix.  A 2x2 `np.ndarray`

    """
    cd = np.array([[hdr[f'CD1_1{key}'], hdr[f'CD1_2{key}']],
                  [hdr[f'CD2_1{key}'], hdr[f'CD2_2{key}']]],
                  dtype=float)
    return cd


def set_cd(hdr, cd, key):
    """
    Function to set a CD matrix to an `astropy.io.fits.Header()`

    Parameters
    ----------
    hdr : `astropy.io.fits.Header()`
       The header to set the data to.

    cd : `np.ndarray()` or list
       The CD matrix as a 2x2 object.

    key : str
       The WCS keyword to set
    """

    hdr[f'CD1_1{key}'] = cd[0, 0]
    hdr[f'CD1_2{key}'] = cd[0, 1]
    hdr[f'CD2_1{key}'] = cd[1, 0]
    hdr[f'CD2_2{key}'] = cd[1, 1]


def get_crval(hdr, key):
    """
    Function to retrieve a CRVAL vector from an `astropy.io.fits.Header()`

    Parameters
    ----------
    hdr : `astropy.io.fits.Header()`
       The header to retrieve from

    key : str
       The WCS key for which WCS to pull from

    Returns
    -------
    crval : `np.ndarray()`
       The CRVAL vector

    """
    crval = np.array([hdr[f'CRVAL1{key}'], hdr[f'CRVAL2{key}']],
                     dtype=float)
    return crval


def set_crval(hdr, crval, key):
    """
    Function to set a CRVAL vector to an `astropy.io.fits.Header()`

    Parameters
    ----------
    hdr : `astropy.io.fits.Header()`
       The header to set the data to.

    crval : `np.ndarray()` or list
       The CRVAL vector as a 2-element object.

    key : str
       The WCS keyword to set
    """

    hdr[f'CRVAL1{key}'] = crval[0]
    hdr[f'CRVAL2{key}'] = crval[1]
