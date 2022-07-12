import numpy as np

from .hdf5table import HDF5Table
from ..utilities import indices


class RDT(HDF5Table):
    COLUMNS=('x','y','lam','val')

    def __init__(self,source,regid,dims=None,**kwargs):
        HDF5Table.__init__(self,dims=dims,**kwargs)
        self.segid=source.segid
        self.regid=regid
        self.pdts={}
        self.pixels=[]

    @property
    def name(self):
        return f'({self.segid},{self.regid})'

    def append(self,pdt):
        ''' append a PDT into this Region Dispersion Table (RDT) '''

        pixel=pdt.pixel
        self.pdts[pixel]=pdt
        self.pixels.append(pixel)

    def decimate(self):
        ''' Decimate (sum over repeated indices) this RDT '''

        if self.pdts and self.dims:
            # extract all the values, but start with any existing data
            x,y=self['x'],self['y']
            lam,val=self['lam'],self['val']
            for pdt in self.pdts.values():
                x.extend(pdt['x'])
                y.extend(pdt['y'])
                lam.extend(pdt['lam'])
                val.extend(pdt['val'])

            # current (aggregated) size of the Table
            if x:
                # ok... clear the space, just to save on memory
                self.pdts.clear()

                # change the data types
                x=np.array(x,dtype=int)
                y=np.array(y,dtype=int)
                lam=np.array(lam,dtype=int)
                val=np.array(val,dtype=float)

                # do the summations
                vv,xx,yy,ll=indices.decimate(val,x,y,lam,dims=self.dims)

                # put these values back in the self
                self.clear()
                self['x']=xx
                self['y']=yy
                self['lam']=ll
                self['val']=vv
