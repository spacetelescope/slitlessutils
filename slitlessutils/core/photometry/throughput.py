import os

import numpy as np
from astropy.io import fits
from scipy.constants import c

from ...config import Config
from ...logger import LOGGER
from .band import Band


class Throughput(Band):

    def __init__(self, *args, **kwargs):
        Band.__init__(self, *args, **kwargs)

        self.zeropoint = kwargs.get('zeropoint', 25.)

    @property
    def photflam(self):
        return self.photfnu * (c / self.photplam) * (1e8 / self.photplam) * 100
        # return 10.**(-0.4*(self.zeropoint+2.408))/self.photplam**2

    @property
    def photfnu(self):
        return 10.**(-0.4 * (self.zeropoint + 48.6))

    @classmethod
    def from_keys(cls, telescope, instrument, band):
        keys = (telescope, instrument, band)

        LOGGER.info(f'loading throughput from keys: {keys}')

        filename = f'{telescope}_{instrument}_{band}.fits'.lower()
        filename = Config().get_reffile(filename, path='bandpasses')

        # check the file is valid.
        if filename is None:
            msg = "Cannot find filter file"
            LOGGER.error(msg)
            raise RuntimeError(msg)

        # read the fits file with the throughput curve
        data, header = fits.getdata(filename, exten=1, header=True)

        obj = cls(data['wavelength'], data['transmission'],
                  unit=header.get('TUNIT1', ''))

        obj.telescope = telescope
        obj.instrument = instrument
        obj.band = band
        obj.filename = filename

        for zlabel in ('ZEROPT', 'MAGZERO', 'ZERO', 'MAG0'):
            if zlabel in header:
                obj.zeropoint = header[zlabel]
                break
        else:
            obj.zeropoint = 25.0
            LOGGER.warning(f"No zeropoint found, using {obj.zeropoint}")

        return obj

    @classmethod
    def from_file(cls, filename):
        tokens = os.path.splitext(filename)
        ext = tokens[-1][1:]
        if ext == 'fits':
            LOGGER.info('reading filter curve from fits file')
            data, header = fits.getdata(filename, exten=1, header=True)
            obj = cls(data['wavelength'], data['transmission'],
                      unit=header.get('TUNIT1', ''))

        elif ext in ('dat', 'txt', 'filt', 'ascii'):
            LOGGER.info("reading filter curve from ascii file")
            wave, tran = np.loadtxt(filename, usecols=(0, 1), unpack=True)
            obj = cls(wave, tran, unit='angstroms')
        else:
            LOGGER.warning(f"File extension {ext} is unknown")
            return

        obj.filename = filename

        return obj

    def __str__(self):
        return f'Throughput curve for {self.telescope}/{self.instrument} {self.band}'
