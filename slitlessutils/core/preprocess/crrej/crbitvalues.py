from ....logger import LOGGER


class BitValues(dict):
    """
    Dictionary to contain the bit values for cosmic rays

    Parameters
    ----------
    None.

    Notes
    -----
    This is just a lookup table
    """

    DEFAULT = 4096

    def __init__(self):
        self['ACS'] = 4096
        self['WFC3'] = 4096

    def __missing__(self, key):
        msg = f"Instrument {key} not found, using default: {self.DEFAULT}"
        LOGGER.warning(msg)
        return self.DEFAULT


BITVALUES = BitValues()
