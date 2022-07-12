from astropy.io import fits
import numpy as np



class SEDFile:
    ''' A class to work with the SED files '''

    def __init__(self,filename,mode='readonly'):
        self.filename=filename
        self.mode=mode

    def __getitem__(self,key):
        if key in self:
            try:
                sed=SED.from_HDU(self._fp[key])
            except:
                filename=self._fp[key].header.get('FILENAME','')
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
