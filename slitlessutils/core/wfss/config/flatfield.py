import numpy as np
from astropy.io import fits

from slitlessutils.logger import LOGGER

from ...utilities import headers


def check_in_field(func):
    """
    A decorator function to check if an (x,y)-coordinate pair is within
    the bounds of a two-dimensional `np.ndarray`.  If the coordinate pair
    is out-of-bounds, the a value of 0.0 is returned.

    Parameters
    ----------
    x : float, int, `np.ndarray`
       The x-coordinates.

    y : float, int, `np.ndarray`
       The y-coordinates.

    l : float, int, `np.ndarray`
       The wavelengths in A.

    Returns
    -------
    f : `np.ndarray`
       The value of the flat field.  In this case, it will always
       be unity.

    """

    def wrapper(self, x, y, l):

        xx = np.rint(x).astype(int)
        yy = np.rint(y).astype(int)

        logic = (xx > 0) & (xx < self.shape[1]) & (yy > 0) & (yy < self.shape[0])
        f = np.where(logic, func(self, xx, yy, l), 0.)
        return f.astype(self.DTYPE)
    return wrapper


class FlatField:
    """
    Base class for all flat field types
    """

    DTYPE = np.float32

    def __init__(self, filename):
        """
        Initializer for flatfield base class

        Parameters
        ----------
        filename : str
            full path to the flat field file
        """
        self.filename = filename

    def update_header(self, hdr):
        """
        Method to update a fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            The fits header to update

        Returns
        -------
        None

        """

        hdr['FFTYPE'] = (self.FTYPE, 'flat field image type')
        hdr['FFNAME'] = (self.filename, 'flat field image name')
        headers.add_stanza(hdr, 'Flat Field Settings', before='FFTYPE')


class UnityFlatField(FlatField):

    """ Class to implement a unity flatfield """

    FTYPE = 'unity'

    def __init__(self):
        FlatField.__init__(self, '')

    def __call__(self, x, y, l):
        """
        Method to call the flat field

        Parameters
        ----------
        x : float, int, `np.ndarray`
           The x-coordinates.

        y : float, int, `np.ndarray`
           The y-coordinates.

        l : float, int, `np.ndarray`
           The wavelengths in A.

        Returns
        -------
        f : `np.ndarray`
           The value of the flat field.  In this case, it will always
           be unity.
        """

        return np.ones_like(x, dtype=self.DTYPE)

    # def update_header(self,hdr):
    #    super().update_header(hdr)


class ImageFlatField(FlatField):
    """
    Class to implement a gray flat field (ie. one with no wavelength
    dependence that likely comes from a direct-image flatfield)
    """

    FTYPE = 'gray'

    def __init__(self, filename):
        FlatField.__init__(self, filename)

    @classmethod
    def from_fits(cls, filename, **kwargs):
        """
        Class Method to load a fits file as a gray flat field

        Parameters
        ----------
        filename : str
            Full path to a valid fits file

        kwargs : dict, optional
            Dictionary of optional keywords passed to the
            `astropy.io.fits.getdata()`

        Returns
        -------
        flat : `FlatField`

        Notes
        -----
        If the input `filename` variable is invalid for any reason, then
        will return a `UnityFlatField`
        """

        obj = cls(filename)

        try:
            if 'ext' not in kwargs:
                kwargs['ext'] = 1
            obj.shape = obj.data.shape
        except BaseException:
            LOGGER.warning("Image flat is invalid, using unity")
            obj = UnityFlatField()

        return obj

    @check_in_field
    def __call__(self, x, y, l):
        """
        Method to call the flat field

        Parameters
        ----------
        x : float, int, `np.ndarray`
           The x-coordinates.

        y : float, int, `np.ndarray`
           The y-coordinates.

        l : float, int, `np.ndarray`
           The wavelengths in A.

        Returns
        -------
        f : `np.ndarray`
           The value of the flat field.  In this case, it will always
           be unity.
        """

        f = self.data[y, x]
        return f

    # def update_header(self,hdr):
    #    super().update_header(hdr)


class PolynomialFlatField(FlatField):
    """
    Class to implement a polynomial flat field of the form:

    .. math::
       f(x,y:\\lambda) = a(x,y) + b(x,y)l + c(x,y)l^2 + d(x,y)l^3+....

    where

    .. math::
       l = \frac{\\lambda - \\lambda_{min}}{\\lambda_{max}-\\lambda_{min}}

    where the coefficients are read from a multi-extension fits file and
    :math:`\\lambda_{min}` and :math:`\\lambda_{max}` are taken from the
    primary header (WMIN and WMAX, respectively).
    """

    FTYPE = 'polynomial'

    def __init__(self, filename):
        FlatField.__init__(self, filename)

    @classmethod
    def from_fits(cls, filename):
        """
        Class Method to load from a fits file

        Parameters
        ----------
        filename : str
            Full path to a valid fits file

        kwargs : dict, optional
            Dictionary of optional keywords passed to the
            `astropy.io.fits.getdata()`

        Returns
        -------
        flat : `FlatField`

        Notes
        -----
        If the input `filename` variable is invalid for any reason, then
        will return a `UnityFlatField`.

        If the fits file has only one valid image extension, then a
        `ImageFlatField` is returned.

        """

        obj = cls(filename)

        with fits.open(obj.filename, mode='readonly') as hdul:
            obj.data = [hdu.data for hdu in hdul if hdu.data is not None]
        obj.order = len(obj.data) - 1

        if obj.order == -1:
            LOGGER.warning("Polynomial flat is invalid, using unity")
            obj = UnityFlatField()
        elif obj.order == 0:
            obj = ImageFlatField.from_fits(filename)
        else:
            obj.wmin = hdul[0].header['WMIN']
            obj.wmax = hdul[0].header['WMAX']
            obj.shape = obj.data[0].shape

        return obj

    @check_in_field
    def __call__(self, x, y, l):
        """
        Method to call the flat field

        Parameters
        ----------
        x : float, int, `np.ndarray`
           The x-coordinates.

        y : float, int, `np.ndarray`
           The y-coordinates.

        l : float, int, `np.ndarray`
           The wavelengths in A.

        Returns
        -------
        f : `np.ndarray`
           The value of the flat field.
        """

        ll = (l - self.wmin) / (self.wmax - self.wmin)
        return sum(f[y, x] * ll**i for i, f in enumerate(self.data))

    # def update_header(self,hdr):
    #    super().update_header(hdr)


def load_flatfield(*args, **kwargs):
    """
    A factory function to load different flat field types

    Parameters
    ----------
    args : tuple
        if no arguments given, then will load a `UnityFlatField`
        if one argument given, then will load either a `ImageFlatField`
           or a `PolynomialFlatField`

    unit : bool, optional
        flag to override args input, and force a `UnityFlatField`

    Returns
    -------
    flat : `FlatField`
        A flatfield object

    Raises
    ------
    NotImplementedError for 2 or more inputs or non-string inputs.

    """

    # sort out the inputs
    unity = kwargs.get('unity', False)

    n = len(args)
    if n == 0 or unity:
        obj = UnityFlatField()
    elif n == 1:
        if isinstance(args[0], str):
            obj = PolynomialFlatField.from_fits(args[0])
        else:
            raise NotImplementedError(f"data type of args: {type(args[1])}")
    else:
        raise NotImplementedError("Unknown number of arguments")

    return obj
