from datetime import datetime

# from ...logger import LOGGER
from ...info import __version__, __code__, __email__, __author__

"""
Methods to work with fits headers
"""


def add_software_log(hdr):
    """
    Method to append to a fits header some properties of the software

    Parameters
    ----------
    hdr : `astropy.io.fits.Header`
        The fits header to update
    """

    hdr.set('CODE', value=__code__, comment='name of software')
    hdr.set('VERSION', value=__version__, comment=f'{__code__} version number')
    hdr.set('AUTHOR', value=__author__, comment=f'author of {__code__}')
    # hdr.set('EMAIL',value=__email__,comment=f'email of author')
    # hdr.set('REF',value=__ref__,
    #        comment=f'publication reference')
    # hdr.set('REFURL',value=__refurl__)

    add_stanza(hdr, "Software Log", before='CODE')


def add_preamble(hdr, **kwargs):
    """
    Method to add a default set of preamble

    Parameters
    ----------
    hdr : `astropy.io.fits.Header`
        The fits header to update

    kwargs : dict, optional
        optional keyword/value pairs to add to the header

    Notes
    -----
    Will add some keywords automatically

    """

    now = datetime.now()
    hdr.set("DATE", value=now.strftime("%Y-%m-%d"),
            comment='date this file was written (yyyy-mm-dd)')
    for k, v in kwargs.items():
        if isinstance(v, tuple):
            n = len(v)
            if n == 1:
                hdr.set(k, value=v[0])
            elif n == 2:
                hdr.set(k, value=v[0], comment=v[1])
            else:
                pass
        else:
            hdr.set(k, value=v)

    add_stanza(hdr, 'File Properties', before='DATE')


def add_stanza(hdr, label, **kwargs):
    """
    Method to add a 'stanza' line to a header

    Parameters
    ----------
    hdr : `astropy.io.fits.Header`
       The fits header to update

    label : str
       The name of the stanza

    kwargs : dict, optional
       Optional parameters passed to `astropy.io.fits.Header().set()`
    """

    hdr.set('', value='', **kwargs)
    hdr.set('', value=f'      / {label}', **kwargs)
    hdr.set('', value='', **kwargs)
