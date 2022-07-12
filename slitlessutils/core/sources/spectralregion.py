import numpy as np

from ..photometry import SED


class SpectralRegion:
    def __init__(self,x,y,w,regid=None):
        self.regid=regid
        self.x=np.array(x,dtype=int)
        self.y=np.array(y,dtype=int)
        self.w=np.array(w,dtype=float)
        self.sed=SED()

    def __len__(self):
        return len(self.x)

    def __iter__(self):
        yield from zip(self.x,self.y,self.w)

    def __getitem__(self,i):
        try:
            value=(self.x[i],self.y[i],self.w[i])
        except:
            value=None
        return value

    def pixels(self):
        yield from zip(self.x,self.y)
