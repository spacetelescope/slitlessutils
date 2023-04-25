from astropy.stats import sigma_clipped_stats
from astropy.io import fits
import os
import numpy as np
from scipy.signal import savgol_filter
import matplotlib.pyplot as plt

from ....config import Config
from ....logger import LOGGER
from ...utilities import headers
from ...wfss import WFSS


def cartesian(sci, src, root, dispaxis=0, inplace=True):
    # if isinstance(data,str):
    #    data=WFSS.observed(data)
    #
    # how to open the fits files
    # if inplace:
    #    mode='update'
    # else:
    #    mode='readonly'

    if dispaxis == 0:
        sci = np.swapaxes(sci, 0, 1)
        src = np.swapaxes(src, 0, 1)
    elif dispaxis == 1:
        pass

    ndisp, ncross = sci.shape
    disp = np.arange(ndisp, dtype=int)
    vals = np.zeros(ndisp, dtype=float)
    stds = np.zeros(ndisp, dtype=float)
    for i in range(ndisp):
        a, m, s = sigma_clipped_stats(sci[i, :], mask=src[i, :])
        vals[i] = m
        stds[i] = s

    smooth = savgol_filter(vals, 33, 3, deriv=0, mode='constant', cval=0.)

    sky = np.tile(smooth, (ncross, 1))

    if dispaxis == 0:
        pass
    elif dispaxis == 1:
        sky = np.swapaxes(sky, 0, 1)

    print(sky.shape, sci.shape)
    fits.writeto(f'{root}.fits', sky, overwrite=True)
    # x=np.arange(ndisp)
    # plt.errorbar(x,vals,stds)
    # plt.plot(x,smooth,color='red')
    # plt.scatter(x,vals)
    # plt.show()
