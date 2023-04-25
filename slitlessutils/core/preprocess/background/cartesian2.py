from astropy.io import fits
import os
import numpy as np


from ....config import Config
from ....logger import LOGGER
from ...utilities import headers
from ...wfss import WFSS


def cartesian(sci, src, dispaxis=0, inplace=True):
    # if isinstance(data,str):
    #    data=WFSS.observed(data)
    #
    # how to open the fits files
    # if inplace:
    #    mode='update'
    # else:
    #    mode='readonly'

    base = os.path.splitext(os.path.basename(data.filename))[0]
    srcfile = f'{base}_src.fits'

    with fits.open(data.filename, mode=mode) as dhdul:
        for hdu in dhdul:
            name = hdu.header.get('EXTNAME')
            if name == 'SCI':

                exten = ('SCI', hdu.header.get('EXTVER', 0))

                sci = dhdul[exten].data
                hdr = dhdul[exten].header

                if os.path.exists(srcfile):
                    src = fits.getdata(srcfile, exten)
                else:
                    src = np.zeros_like(src, dtype=int)

                if dispaxis == 0:
                    pass
                elif dispaxis == 1:
                    sci = np.swapaxes(sci, 0, 1)
                    src = np.swapaxes(src, 0, 1)

                ndisp, ncross = sci.shape
                disp = np.arange(ndisp, dtype=int)
                vals = np.zeros(ndisp, dtype=float)

                for i in range(ndisp):
                    a, m, s = sigma_clipped_stats(sci[i, :], mask=sky[i, :])
                    vals[i] = m

                smooth = savgol_filter(vals, 37, 3, deriv=0, mode='constant', cval=0.)

                sky = np.tile(smooth, (ncross, 1))
                print(sky.shape, sci.shape)
