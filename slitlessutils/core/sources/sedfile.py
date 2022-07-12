from astropy.io import fits
import numpy as np

from ...logger import LOGGER
from ..photometry import SED

class SEDFile:
    ''' A class to work with the SED files '''

    def __init__(self,filename,mode='readonly'):
        self.filename=filename
        self.mode=mode

    def __getitem__(self,key):
        if isinstance(key,(tuple,list)) and len(key)==2:
            exten=(str(key[0]),int(key[1]))
        elif isinstance(key,(int,np.integer)):
            exten=(str(key),0)
        elif isinstance(key,str):
            exten=(key,0)
        else:
            exten=None
            
        if exten:
            try:
                sed=SED.from_HDU(self._fp[exten])
            except:
                filename=self._fp[exten].header.get('FILENAME','')
                sed=SED.from_file(filename)
        else:
            sed=SED()
        return sed

    def __contains__(self,k):
        return k in self._fp

    def __enter__(self):
        self._fp=fits.open(self.filename,mode=self.mode)
        return self

    def __exit__(self,etype,eval,etrace):
        self.close()

    def close(self):
        self._fp.close()
