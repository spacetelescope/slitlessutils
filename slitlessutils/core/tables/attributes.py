import numpy as np

"""
functions to load and write attributes to an HDF5 file.  Attributes
are meta data that are conceptually similar to the fits header keywords
"""


def load(h5, key, ptype=None):
    """
    Function to load attributes from the HDF5 file

    Parameters
    ----------
    h5 : a valid `h5py` object (can be dataset or group)
        The object to load attributes from

    key : str
        The name of the keyword to use

    ptype : type or None, optional
        The typing will be done automatically, but this is a way of
        overriding it.  If None, then use the automatic typing.  Default
        is None

        The rules for automatic typing are:
        1) if the value is byte data then:
           decode as UTF-8 and make all lowercase
             a) if 'true', then return True
             b) if 'false', then return False
             c) if 'none', then return None
             else) return the UTF-8 decoded result
        2) if value is a float then
           a) if value is nan, then return None
           else) return the float value
        3) if none of these, then return the value as is

        The ptype function applied as it was being returned

    Return
    ------
    val : arb. type
        The returned typed data

    """

    val = h5.attrs.get(key, None)
    if isinstance(val, bytes):
        val = val.decode('UTF-8')
        vlo = val.lower()
        if vlo == 'true':
            val = True
        elif vlo == 'false':
            val = False
        elif vlo == 'none':
            val = None
        else:
            pass
    elif isinstance(val, float):
        if np.isnan(val):
            val = None
        else:
            pass
    else:
        pass

    if isinstance(ptype, type):
        val = ptype(val)

    return val


def write(h5, key, val):
    """
    Function to write an attribute to an HDF5 object

    Parameters
    ----------
    h5 : a valid `h5py` object (dataset or group)
        The HDF5 object to write an attribute to

    key : str
        The name of the keyword to write

    val : arbitrary type
        The data to put into the attribute.

    """

    if val is not None:
        if isinstance(val, (bool, str)):
            h5.attrs[key] = np.bytes_(val)
        else:
            h5.attrs[key] = val
