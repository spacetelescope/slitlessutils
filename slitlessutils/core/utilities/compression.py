# import subprocess
import os
import gzip
import shutil

from ...logger import LOGGER

"""
Methods to work with gzip utiltiles
"""

SUPPORTED_TYPES = ('gz',)


def compress(filename, comptype='gz', keeporig=False):
    """
    Method to compress a file

    Parameters
    ----------
    filename : str
        Name of the file to compress

    comptype : str, optional
        Type of the compression.  Default is 'gz'

    keeporig : bool, optional
        Flag to keep the original file as well.  Default is False

    Returns
    -------
    newfile : str
        name of the zipped file
    """

    if comptype not in SUPPORTED_TYPES:
        LOGGER.warning(f'Unsupported compression type: {comptype}')
        return None

    zipfile = filename+'.'+comptype

    with open(filename, 'rb') as new_file:
        if comptype in ('gz', 'gzip'):
            with gzip.open(zipfile, 'wb') as zip_file:
                zip_file.writelines(new_file)

    if not keeporig:
        os.remove(filename)

    return zipfile


def uncompress(filename, keeporig=False):
    """
    Method to uncompress a file

    Parameters
    ----------
    filename : str
        Name of the file to uncompress

    keeporig : bool, optional
        Flag to keep the original file.  Default is False

    Returns
    -------
    newfile : str
        uncompressed name of the file
    """

    newfile, ext = os.path.splitext(filename)
    comptype = ext[1:]

    if comptype not in SUPPORTED_TYPES:
        LOGGER.warning(f'Unsupported compression type: {comptype}')
        return

    with open(newfile, 'wb') as orig_file:
        if comptype in ('gz', 'gzip'):
            with gzip.open(filename, 'rb') as zip_file:
                shutil.copyfileobj(zip_file, orig_file)

    if not keeporig:
        os.remove(filename)

    return newfile


if __name__ == '__main__':
    with open('t.txt', 'w') as f:
        print('this is a file', file=f)
    compress('icoi3qcdq_flt.fits')
    uncompress('icoi3qcdq_flt.fits.gz')
