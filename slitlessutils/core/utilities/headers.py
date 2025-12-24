"""
Method to work with fits headers
"""

from datetime import datetime

from .get_metadata import get_metadata


def add_software_log(hdr):
    """
    Method to append to a fits header some properties of the software

    Parameters
    ----------
    hdr : `astropy.io.fits.Header`
        The fits header to update
    """

    # get the package meta data
    meta = get_metadata()
    package = meta.get('Name', '')

    # update the header
    hdr.set('PACKAGE', value=package, comment='software package')

    hdr.set('VERSION', value=meta.get('Version', ''),
            comment=f'{package} version number')

    hdr.set('AUTHOR', value=meta.get('Author', ''),
            comment=f'{package} author')

    # finalize by adding a stanza header
    add_stanza(hdr, "Software Log", before='PACKAGE')


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
