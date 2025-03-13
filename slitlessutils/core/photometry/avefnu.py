import numpy as np
from astropy.utils import minversion

from ...logger import LOGGER


def avefnu(sed, band):
    """
    Function to compute the bandpass-averaged flux of a spectrum given by:

    .. math::
       \\left<F_\nu\right> = \frac{\\int F_\nu T_\nu d\nu/\nu}{\\int T_\nu d\nu/\nu}

    where :math:`F_\nu` is the spectrum, :math:`T_{\nu}` is the transmission
    curve and :math:`\nu` is frequency.

    Parameters
    ----------
    sed : `slitlessutils.core.photometry.SED`
        An `SED` object

    band : `slitlessutils.core.photometry.Band`
        A `Band` object

    Returns
    -------
    ave : float
        The bandpass-averaged flux

    Notes
    -----
    Uses a trapezoidal rule for numerical integration, which is fine for
    well sampled SEDs and Bandpasses.  This could get problematic for poorly
    sampled things

    """

    ave = np.nan
    if sed.wmin > band.wmax:
        LOGGER.warning(f"Bandpass {band.name} is entirely too blue")
    elif band.wmin > sed.wmax:
        LOGGER.warning(f"Bandpass {band.name} is entirely too red")
    elif (sed.wmin > band.wmin) or (sed.wmax < band.wmax):
        LOGGER.warning(f"Bandpass {band.name} does not cover full range of sed")
    else:
        # import matplotlib.pyplot as plt
        # plt.plot(band.wave,band.tran*np.amax(sed.flam)/np.amax(band.tran))
        # plt.plot(sed.lamb,sed.flam)
        # plt.xlim(1200,2000)
        # plt.show()

        fnu = sed(band.wave, fnu=True)

        if minversion(np, '2.0'):
            trap_function = np.trapezoid
        else:
            trap_function = np.trapz
        ave = trap_function(fnu * band.tran / band.freq,
                            x=band.wave) / band.fnunorm

    return ave
