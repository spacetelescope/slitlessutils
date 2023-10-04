from .affinetweak import AffineTweak
from ....logger import LOGGER


def upgrade_wcs(reffile, inputs, key0='A', key1='', newfile=None, inplace=False):
    """
    Method to upgrade the WCS in one fits image to the equivalent in another
    using an affine transformation.  For more information on this, see the
    class ``AffineTweak`` for more information.  This is usually necessary for
    HST data, as at some point the WCS in the direct image was updated with
    respect to Gaia, but the WFSS images were not.

    Parameters
    ----------
    reffile : str
        Full path to the file that the tweak is computed from.

    inputs : str or list or tuple
        If this is a str, then it is the name of the file to be upgrade.  If
        it is a list or tuple, then it is many such filenames.

    key1 : str, optional
        The WCS key for the final astrometry.  Default is ''.

    key0 : str, optional
        The WCS key for the final astrometry.  Default is 'A'.

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

    # compute the tweak
    tweak = AffineTweak(reffile, key0=key0, key1=key1)

    if isinstance(inputs, str):
        return tweak(inputs, inplace=inplace, newfile=newfile)
    elif isinstance(inputs, (list, tuple)):
        return [tweak(i, inplace=inplace, newfile=newfile) for i in inputs]
    else:
        LOGGER.warning(f"Invalid datatype ({type(inputs)}) for input.")
        return
