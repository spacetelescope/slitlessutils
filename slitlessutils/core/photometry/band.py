import numpy as np
from scipy.constants import c

from ...logger import LOGGER

class Band:
    ''' A base class to handle transmission curves and sensitivity files '''

    def __init__(self,wave,tran,where=None,unit='angstrom'):
        if where is None:
            where = np.arange(len(wave),dtype=int)

        self.unit=unit.lower()
        self.wave=wave
        self.tran=tran

        if self.unit in ('micron','microns','um'):
            self.wave*=1e4
        elif self.unit in ('angstrom','angstroms','a',''):
            pass
        else:
            LOGGER.warning(f'Unknown wavelength units {self.unit}. Assuming A')

        self.freq=(c/self.wave)*1e10
        
        # set ranges
        self.wmin=np.amin(self.wave[where])
        self.wmax=np.amax(self.wave[where])

        # compute the pivot wavelength
        num=np.trapz(self.tran[where],x=self.wave[where])
        den=np.trapz(self.tran[where]/self.wave[where]**2,x=self.wave[where])
        self.photplam=np.sqrt(num/den)

        # compute the max value
        self.tmax=np.amax(self.tran)

        # compute normalization
        self.fnunorm=np.trapz(self.tran/self.freq,x=self.freq)

        
        
    def __call__(self,l,left=0.,right=0.):
        s=np.interp(l,self.wave,self.tran,left=left,right=right)
        return s

    def __mul__(self,a):
        self.tran*=a
        self.tmax*=a
        return self

    def __rmul__(self,a):
        return self.__mul__(a)

    def __imul__(self,a):
        return self.__mul__(a)
