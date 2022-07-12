import numpy as np
from scipy.constants import c


from ...logger import LOGGER
def avefnu(sed,band):
    
    if sed.wmin > band.wmax:
        LOGGER.warning(f"Bandpass {band.name} is entirely too blue")
    elif band.wmin > sed.wmax:
        LOGGER.warning(f"Bandpass {band.name} is entirely too red")
    elif (sed.wmin > band.wmin) or (sed.wmax < band.wmax):
        LOGGER.warning(f"Bandpass {band.name} does not cover full range of sed")



    fnu=sed(band.wave,fnu=True)
    ave=np.trapz(fnu*band.tran/band.freq,x=band.freq)/band.fnunorm
    return ave
    
