from astropy.io import fits
import numpy as np

from ..photometry import SED


class SEDFile:
    """
    A class to load SEDs used for simulating WFSS images.

    Instatiates a variable that operates like a dictionary, whose keywords
    refer to the ``EXTNAME`` and ``EXTVER``, which denote the ``SEGID``
    and ``REGID``, respectively.  Therefore there are two possible
    options for the keyword:

    (1) scalar value : this is interpreted as the ``EXTNAME`` keyword, but
        must be a datatype that can cast as a `str`.

    (2) list or tuple : this is required to have only two elements, which
        refer to the ``EXTNAME`` and ``EXTVER``, respectively.

    In all cases, this will return a ``~slitlessutils.photometry.SED``
    object.


    Parameters
    ----------
    filename : str
        A string variable that indicates the name of a multi-extension fits
        file, where each extension is a `~astropy.io.fits.BinTableHDU`
        that contains a fits table of `~astropy.table.Table`.

    mode : str, optional
        A string with the same meaning as in `~astropy.io.fits.open`.

    Notes
    -----
    This will operate like a dict, but does not have a dict internally.

    """

    def __init__(self, filename, mode='readonly'):
        self.filename = filename
        self.mode = mode

    def __getitem__(self, key):
        """
        Method to read an SED from the file

        Parameters
        ----------
        key : many possible types
           Key used to extract a spectrum from the file.  The rules for
           the possible types of key are:

           if key is a list or tuple and has 2 elements, then key is
               interpreted as the EXTNAME and EXTVER for reading from file.

           if key is int or str then it is interpreted as the EXTNAME

           All other types will return an empty `SED`


        Returns
        -------
        sed : `su.core.photometry.SED`
           The SED
        """

        if isinstance(key, (tuple, list)) and len(key) == 2:
            exten = (str(key[0]), int(key[1]))
        elif isinstance(key, (int, np.integer)):
            exten = (str(key), 0)
        elif isinstance(key, str):
            exten = (key, 0)
        else:
            exten = None

        if exten:
            try:
                sed = SED.from_HDU(self._fp[exten])
            except BaseException:
                filename = self._fp[exten].header.get('FILENAME', '')
                sed = SED.from_file(filename)
        else:
            sed = SED()
        return sed

    def __contains__(self, k):
        """
        Method to implement the 'in' operator
        """

        return k in self._fp

    def __enter__(self):
        """
        Method to implement the context manager
        """

        self._fp = fits.open(self.filename, mode=self.mode)
        return self

    def __exit__(self, etype, eval, etrace):
        """
        Method to implement the context manager
        """
        self.close()

    def close(self):
        """
        Method to explicitly close the file.
        """
        self._fp.close()
