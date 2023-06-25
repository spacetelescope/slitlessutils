from astropy.io import fits
import numpy as np
from scipy.constants import c
import socket


import os
# from urllib import request
import requests

from .throughput import Throughput
from .avefnu import avefnu
from ...logger import LOGGER
from ...config import Config


class SED:
    """
    Class for SED functionalities.  Emulates list behaviors.

    Parameters
    ----------
    args : tuple
        See the class method `SED.append()` for more details

    n : int, optional
        The initial size for an empty SED.  This is to govern the dynamic
        state of the internal numpy array.  Default is 100

    kwargs : dict, tuple
        See the class method `SED.append()` for more details

    Notes
    -----
    1) Internally, the data are stored in a structured array whose size
       is larger than the data it holds, so if data are appended, then
       the array does not need to change.  If the array's size does need
       to increase, then there are parameters to control how fast it
       grows.

    2) Alternatively, data structures like pandas (as they claim poor
       performance on small arrays and/or appending) or lists (as the data
       are needed as `np.ndarray` for calculations).  Therefore, this
       approach was used to facilitate real-time appending and data
       accessing.  But this does require several dunder-methods to emulate
       the dictionary-like nature.

    3) The SED object can be accessed either as a dictionary or by attribute.

    """

    # dtypes for a structured array
    DTYPE = [('lamb', np.float32),
             ('flam', np.float32),
             ('func', np.float32),
             ('cont', np.float32),
             ('npix', np.uint16)]
    FORMAT = ('%.4e', '%.4e', '%.4e', '%.4e', '%6u')
    GROW = 2     # factor to resize a dynamic array

    def __init__(self, *args, n=100, **kwargs):
        self.reset(*args, n=n, **kwargs)

    def __len__(self):
        """
        Method to override the len() function
        """
        return self.count

    def __getitem__(self, key):
        """
        Method to retrieve a column

        Parameters
        ----------
        key : str
            The column name to retrieve

        Results
        -------
        v : `np.ndarray`
            The column data
        """
        return self._data[key][:self.count]

    def __setitem__(self, key, value):
        """
        Method to set data to a column

        Parameters
        ----------
        key : str
            The column name to set

        value : tuple, list, or `np.ndarray`
            The data to set to the column

        """
        self._data[key][:self.count] = value

    @property
    def data(self):
        """
        Property to contain the data as a structured array

        Returns
        -------
        data : `np.ndarray`
           The structured array
        """
        return self._data[:self.count]

    @data.setter
    def data(self, d):
        """
        Setter to set data

        Parameters
        ----------
        d : `np.ndarray`
            The data to set
        """
        self._data = d

    # def normalize(self,band,fnu):
    def normalize(self, wave, flux, abmag=False):
        """
        Method to renormalize the spectrum

        Parameters
        ----------
        wave : int,float, or `su.core.photometry.Throughput`
           This is a flexible object that dictates where the spectrum should
           be normalized.  If it is passed as a numeric type, then it will
           be used to interpolate the spectrum *AT* that wavelength.  If it
           is passed as a `su.core.photometry.Throughput`, then it will be
           used to normalize through a filter.

        flux : int, float
           The value to force the spectrum to have.  This can be either
           fnu-based flux or an AB magnitude, as controled by the `abmag`
           optional parameter.

        abmag : bool, optional
           A flag that determines how the `flux` parameter is to be
           interpreted.  If `abmag==True`, then consider flux as an
           AB magnitude.  If `abmag==False`, then consider flux as a flux
           in fnu units.

        """

        if isinstance(wave, (float, int)):
            den = self(wave)
        elif isinstance(wave, Throughput):
            den = avefnu(self, wave)
        else:
            LOGGER.warning(f'Type of wave is unsupported {type(wave)}.')
            return

        if abmag:
            num = 10.0**(-0.4*(flux+48.6))
        else:
            num = flux

        self *= (num/den)

    def reset(self, *args, n=100, **kwargs):
        """
        Method to reset the state of SED object

        Parameters
        ----------
        args : tuple
           See the class method `SED.append()` for more details

        n : int, optional
           The initial size for an empty SED.  This is to govern the dynamic
           state of the internal numpy array.  Default is 100

        kwargs : dict, tuple
           See the class method `SED.append()` for more details

        """

        self.count = 0
        if len(args) == 2 and np.shape(args[0]) == np.shape(args[1]):
            self.data = np.zeros(self.GROW*len(args[0]), dtype=self.DTYPE)
            self.append(*args, **kwargs)
        else:
            self.data = np.zeros(n, dtype=self.DTYPE)

    def append(self, lamb, flam, func=None, cont=None, npix=None):
        """
        Method to append new elements to this SED

        Parameters
        ----------
        lamb : int, float, `np.ndarray`
           The wavelength elements.

        flam : int, float, `np.ndarray`
           The flux elements, in flam units

        func : int, float, `np.ndarray` or `None`, optional
           The uncertainty elements in flam units.  Default is None

        cont : int, float, `np.ndarray` or `None`, optional
           The contamination estimate.  Default is None.

        npix : int, `np.ndarray`, or `None`, optional
           The number of pixels in this coordinate.  Default is None

        Notes
        -----
        1) The shape of `lamb` and `flam` should match.
        2) If `func`, `cont`, and/or `npix` are not None, then their
           shape should match lamb. If they are specified and do not match
           the shape, then a NaN or zero is used, as applicable.

        """

        lshape = np.shape(lamb)
        if np.shape(flam) != lshape:
            LOGGER.warning('Shape of lamb and flam do not match.')
            return
        else:
            if lshape:
                nadd = lshape[0]
            else:
                nadd = 1
        length = len(self._data)

        # check to grow the array
        nneed = nadd-length-self.count
        if nneed > 0:    # grow the array
            data = np.zeros(length+self.GROW*nneed, dtype=self.DTYPE)
            data[:self.count] = np.copy(self.data)
            self._data = np.copy(data)

        # check the inputs
        if func is None or np.shape(func) != lshape:
            func = np.full_like(lamb, np.nan)
        if cont is None or np.shape(cont) != lshape:
            cont = np.full_like(lamb, np.nan)
        if npix is None or np.shape(npix) != lshape:
            npix = np.zeros_like(lamb, dtype=int)

        # put data in array
        self._data['lamb'][self.count:self.count+nadd] = lamb
        self._data['flam'][self.count:self.count+nadd] = flam
        self._data['func'][self.count:self.count+nadd] = func
        self._data['cont'][self.count:self.count+nadd] = cont
        self._data['npix'][self.count:self.count+nadd] = npix
        self.count += nadd

    def __bool__(self):
        """
        Method to enable to boolean testing based on the count
        """

        return self.count > 0

    def __mul__(self, a):
        """
        Enable multiplication, scaling the spectrum by a constant

        Parameters
        ----------
        a : float, int
           Scalar constant

        Returns
        -------
        sed : `SED`
           A scaled spectrum


        """

        sed = SED()
        sed.count = self.count
        sed._data = self._data.copy()

        sed['flam'] *= a
        sed['func'] *= a
        sed['cont'] *= a
        return sed

    def __rmul__(self, a):
        """
        Enable multiplication, scaling the spectrum by a constant

        Parameters
        ----------
        a : float, int
           Scalar constant

        Returns
        -------
        sed : `SED`
           A scaled spectrum

        """
        return self.__mul__(a)

    def __imul__(self, a):
        """
        Enable self-multiplication, scaling the spectrum by a constant

        Parameters
        ----------
        a : float, int
           Scalar constant

        Returns
        -------
        sed : `SED`
           A scaled version of the self spectrum

        """
        if np.isinf(a):
            raise ValueError("Cannot multiply by infinity")

        self['flam'] *= a
        self['func'] *= a
        self['cont'] *= a
        return self

    def __add__(self, a):
        """
        Enable addition, shifting the spectrum by a constant

        Parameters
        ----------
        a : float, int, or `SED`
           Scalar constant

        Returns
        -------
        sed : `SED`
           The shifted spectrum
        """
        sed = SED()
        if isinstance(a, SED):
            if np.allclose(a.lamb, self.lamb):
                sed.count = self.count
                sed._data = self._data.copy()
                sed.flam += a.flam
        else:
            sed.flam += a

        return sed

    def __radd__(self, a):
        """
        Enable addition, shifting the spectrum by a constant

        Parameters
        ----------
        a : float, int, or `SED`
           Scalar constant

        Returns
        -------
        sed : `SED`
           The shifted spectrum
        """
        return self.__add__(a)

    def __iadd__(self, a):
        """
        Enable addition, shifting the spectrum by a constant

        Parameters
        ----------
        a : float, int, or `SED`
           Scalar constant

        Returns
        -------
        sed : `SED`
           The shifted self-spectrum
        """

        if isinstance(a, SED):
            if np.allclose(a.lamb, self.lamb):
                self.flam += a.flam
            else:
                raise RuntimeError
        else:
            self.flam += a
        return self

    def redshift(self, z):
        """
        Method to redshift a spectrum

        Parameters
        ----------
        z : float, int
           The redshift

        """

        self['lamb'] *= (1+z)

    @property
    def wmin(self):
        """
        Property for the wavelength minimum
        """

        return np.amin(self.lamb)

    @property
    def wmax(self):
        """
        Property for the wavelength maximum
        """

        return np.amax(self.lamb)

    def __getattr__(self, k):
        """
        Method to enable attribute-like access.

        Parameters
        ----------
        k : str
           The attribute to get

        Returns
        -------
        v : any type
           The attribute value

        Notes
        -----
        If the keyword is one of the data structures keys, then return
        that element of the data.  Otherwise, call __getattr__ like normal.
        """

        if k in ('lamb', 'flam', 'func', 'cont', 'npix'):
            return self[k]
        else:
            return super().__getattr__(k)

    def __setattr__(self, k, v):
        """
        Method to enable attribute-like setting

        Parameters
        ----------
        k : str
           The attribute to get

        v : any type
           The attribute value

        Notes
        -----
        If the keyword is one of the data structures keys, then return
        that element of the data.  Otherwise, call __setattr__ like normal.
        """

        if k in ('lamb', 'flam', 'func', 'cont', 'npix'):
            self[k] = v
        else:
            super().__setattr__(k, v)

    def __call__(self, wave, fnu=False, **kwargs):
        """
        Method to evaluate the spectrum via interpolation

        Parameters
        ----------
        wave : float, int, or `np.ndarray`
           The wavelength values

        fnu : bool, optional
           A flag to convert the flam spectrum to fnu.  Default is False

        kwargs : dict, optional
           Optional keywords passed to `np.interp()`

        Returns
        -------
        flux : float or `np.ndarray`
           The spectral values at the input wavelengths
        """

        if not self:
            return np.zeros_like(wave, dtype=float)

        g = np.where(np.isfinite(self.flam))[0]
        lamb = self.lamb[g]
        flam = self.flam[g]
        flux = np.interp(wave, lamb, flam, **kwargs, left=flam[0], right=flam[-1])

        # flux=np.interp(wave,self.lamb[g],self.flam[g],**kwargs,
        #               left=self.flam[g[0]],right=self.flam[g[-1]])

        if fnu:
            flux *= ((wave/c)*(wave/1e10))

        return flux

    @classmethod
    def from_HDU(cls, hdu):
        """
        Classmethod to load a spectrum from a header-data unit (HDU)

        Parameters
        ----------
        hdu : `astropy.io.fits.HDU`
           The HDU to get the data from

        Returns
        -------
        sed : `SED`
           The SED
        """

        # parse data from the HDU
        if 'uncertainty' in hdu.data.names:
            func = hdu.data['uncertainty']
        else:
            func = None

        try:
            lamb, flux = hdu.data['wavelength'], hdu.data['flux']
        except BaseException:
            lamb, flux = hdu.data['lamb'], hdu.data['flam']

        # adjust units of wavelengths
        units = hdu.header.get('TUNIT1', '').lower()
        if units in ('angstrom', 'angstroms', 'a', ''):
            pass
        elif units in ('micron', 'microns', 'um'):
            lamb *= 1e4
        else:
            LOGGER.warning(f'Unknown wavelength units {units}. Assuming A')

        # adjust units of flux density
        units = hdu.header.get("TUNIT2", '').lower()
        if units in ('fnu',):
            const = (c/lamb)*(1e10/lamb)
            flam = flux*const
            flamunc = func*const
        elif units in ('flam', 'erg/(s*cm**2*aa)'):
            flam = np.copy(flux)
            flamunc = np.copy(func)
        else:
            flam = np.copy(flux)
            flamunc = np.copy(func)
            LOGGER.warning(f'Unknown flux units {units}.  Assuming Flam')

        obj = cls(lamb, flam, func=flamunc)

        # check some things in the header
        if 'REDSHIFT' in hdu.header:
            obj.redshift(hdu.header['REDSHIFT'])
        return obj

    @classmethod
    def from_file(cls, filename, exten=1):
        """
        Classmethod to load an SED from a file

        Parameters
        ----------
        filename : str
           The name of the file to load

        exten : int or tuple or string, optional
           The extension in if passed a fits file. Default is 1

        Returns
        -------
        sed : SED
           The SED object

        Notes
        -----
        The file type is interpreted from the extension on the filename
        (ie. the token following the last '.').  Valid file types are:

        |-----------------|-------------------------------------------|
        | File extension  | Interpretation                            |
        |-----------------|-------------------------------------------|
        | dat txt ascii   | space-delimited columnar data of the form |
        | sed or blank    | lamb,flam                                 |
        |                 |                                           |
        |-----------------|-------------------------------------------|
        | fits, fit       | a Binary Table fits file. The columns are |
        |                 | specified in `from_HDU()`                 |
        |-----------------|-------------------------------------------|
        | csv             | a comma-separated value file, where the   |
        |                 | first row specifies the columns, valid    |
        |                 | columns are: 'lamb','flam', 'uncertainty' |
        |                 | 'contamination', and 'npixel'             |
        |-----------------|-------------------------------------------|

        """

        if os.path.exists(filename):
            ext = os.path.splitext(filename)[1][1:].lower()

            if ext in ('dat', 'txt', 'ascii', 'sed', ''):
                lamb, flam = np.loadtxt(filename, usecols=(0, 1), unpack=True)
                obj = cls(lamb, flam)

            elif ext in ('fits', 'fit'):
                # exten=kwargs.get('exten',1)
                with fits.open(filename, mode='readonly') as hdul:
                    obj = cls.from_HDU(hdul[exten])
                return obj

            elif ext in ('csv',):
                raise NotImplementedError('csv file not supported')

            else:
                LOGGER.warning(f"File type {ext} is not supported")
                return None
        else:
            LOGGER.warning(f"SED File not found: {filename}")
            return None

        obj.filename = filename      # record the filename
        return obj

    @staticmethod
    def get_from_CDBS(atlas, filename):
        """
        Staticmethod to retrieve a spectrum from the Calibration Database
        System (CDBS)

        Parameters
        ----------
        atlas : str
            The spectral atlas.

        filename : str
            The filename in the atlas

        """

        # first check that internet is alive
        IP = socket.gethostbyname(socket.gethostname())
        if IP == '127.0.0.1':
            LOGGER.warning(f'Invalid IP: {IP}.  Skipping download from CDBS.')
            return

        # base URL for CDBS
        url = 'https://archive.stsci.edu/hlsps/reference-atlases/cdbs/grid/'

        # do something for each valid atlas
        if atlas == 'bc95':
            serverfile = f'{url}{atlas}/templates/{filename}'
        elif atlas == 'pickles':
            serverfile = f'{url}{atlas}/dat_uvk/{filename}'
        else:
            LOGGER.warning(f'Invalid CDBS Atlas: {atlas}.  Skipping download.')
            return

        # try getting the file in a secure way
        resp = requests.get(serverfile, timeout=5)
        if resp == 404:
            LOGGER.warning(f'Cannot find server-side file: {serverfile}')
            return
        else:
            try:
                with open(filename, 'wb') as fp:
                    fp.write(resp.content)
                return filename
            except BaseException:
                LOGGER.warning(f"Cannot write local file: {filename}")

        # try getting the file
        # try:
        #     f, r = request.urlretrieve(serverfile, filename)
        #     localfile = filename
        #     return localfile
        # except BaseException:
        #     LOGGER.warning(f'Cannot find server-side file: {serverfile}')
        #     return

    @classmethod
    def from_CDBS(obj, atlas, filename, cleanup=True):
        """
        Classmethod to load an `SED` from a spectrum in the Calibration
        Database System (CDBS)

        Parameters
        ----------
        atlas : str
            The spectral atlas in CDBS

        filename : str
            The filename in the atlas

        cleanup : bool, optional
            Flag to delete the local copy of the CDBS file when finished.
            Default is True

        Results
        -------
        sed : `SED`
            The SED.
        """

        obj.get_from_CDBS(atlas, filename)

        sed = obj.from_file(filename)
        if cleanup:
            os.remove(filename)

        return sed

    def write_file(self, filename, **kwargs):
        """
        Method to write the SED to a file

        Parameters
        ----------
        filename : str
            The name of the file to write.  The file format is taken from the
            extension specified in this string.

        kwargs : dict, optional
            See the `as_HDU()` for these keywords

        Notes
        -----
        The file format is determined from the passed filename's extension:

        |-----------------|-------------------------------------------|
        | File extension  | Interpretation                            |
        |-----------------|-------------------------------------------|
        | dat txt ascii   | space-delimited columnar data of the form |
        | sed, or blank   | lamb,flam                                 |
        |                 |                                           |
        |-----------------|-------------------------------------------|
        | fits, fit       | a Binary Table fits file. The columns are |
        |                 | specified in `as_HDU()`                   |
        |-----------------|-------------------------------------------|
        | csv             | a comma-separated value file, where the   |
        |                 | first row specifies the columns, valid    |
        |                 | columns are: 'lamb','flam', 'uncertainty' |
        |                 | 'contamination', and 'npixel'             |
        |-----------------|-------------------------------------------|

        """

        # check if path exists
        path = os.path.dirname(filename)
        if not os.path.isdir(path):
            os.mkdir(path)

        ext = os.path.splitext(filename)[1][1:].lower()
        if ext in ('csv',):
            delimiter = ','
            np.savetxt(filename, self.data, delimiter=delimiter, fmt=self.FORMAT,
                       header=delimiter.join(self._data.dtype.names))
        elif ext in ('txt', 'dat', 'sed', 'ascii', ''):
            delimiter = ' '
            np.savetxt(filename, self.data, delimiter=delimiter, fmt=self.FORMAT,
                       header=delimiter.join(self._data.dtype.names))
        elif ext in ('fits', 'fit'):
            hdu = self.as_HDU(**kwargs)
            hdu.writeto(filename, overwrite=True)
        else:
            delimiter = ' '
            np.savetxt(filename, self.data, delimiter=delimiter, fmt=self.FORMAT,
                       header=delimiter.join(self._data.dtype.names))

    def as_HDU(self, **kwargs):
        """
        Method to package the spectrum in to an `astropy.io.fits.BinTableHDU`

        Parameters
        ----------
        kwargs : dict, optional
            Keyword/value pairs set as header keyword/value pairs

        """

        fluxunits = Config().fluxunits

        hdu = fits.BinTableHDU(self.data)

        hdu.header.set('TUNIT1', value='angstrom', after='TFORM1')
        hdu.header.set('TUNIT2', value=fluxunits, after='TFORM2')
        hdu.header.set('TUNIT3', value=fluxunits, after='TFORM3')
        hdu.header.set('TUNIT4', value=fluxunits, after='TFORM4')
        hdu.header.set('TUNIT5', value='number', after='TFORM5')

        for k, v in kwargs.items():
            hdu.header[k] = v
        return hdu


if __name__ == '__main__':
    x = SED()
    l = np.arange(50)
    f = np.arange(1, 51)
    x.append(l, f)
    y = x*2.
    print(y.data)
