import subprocess
import os

from ...logger import LOGGER

"""
Methods to work with gzip utiltiles
"""


def gzip(filename):
    """
    Method to gzip a file

    Parameters
    ----------
    filename : str
        Name of the file to gzip

    Returns
    -------
    newfile : str
        Gzipped file name
    """

    r = subprocess.call(['gzip', '-f', filename])
    if r != 0:
        LOGGER.warning(f'Cannot gzip {filename} with error: {r}')
    return filename+'.gz'


def gunzip(filename):
    """
    Method to gunzip a file

    Parameters
    ----------
    filename : str
        Name of the file to gunzip

    Returns
    -------
    newfile : str
        Gunzipped file name
    """

    r = subprocess.call(['gunzip', filename])
    if r != 0:
        LOGGER.warning(f'Cannot gunzip {filename} with error: {r}')

    newfilename, exten = os.path.splitext(filename)

    return newfilename
