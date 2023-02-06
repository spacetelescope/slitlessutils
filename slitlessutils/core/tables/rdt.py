import numpy as np

from .hdf5table import HDF5Table
from ..utilities import indices


class RDT(HDF5Table):
    """
    Class to hold the region-dispersion table (RDT).  Here the word `region`
    refers to a spectral region contained *inside* a source.  This
    is largely for extended sources.

    inherits from `HDF5Table`

    """


    # the columns for this table    
    COLUMNS=('x','y','lam','val')
    
    def __init__(self,source,regid,dims=None,**kwargs):
        """
        Initializer
        
        Parameters
        ----------
        source : `su.core.sources.Source`
            The source to load

        regid : int
            The region ID (recalling that the region is a spectral region 
            contained within the source). 

        dims : tuple or None, optional
            Value passed to `HDF5Table`

        kwargs : dict, optional
            Key value pairs to pass to `HDF5Table`
        """


        HDF5Table.__init__(self,dims=dims,**kwargs)
        self.segid=source.segid
        self.regid=regid
        self.pdts={}
        self.pixels=[]

    @property
    def name(self):
        return f'({self.segid},{self.regid})'

    def append(self,pdt):
        """
        Method to append another pdt into this table

        Parameters
        ----------
        pdt : `su.core.tables.PDT`
            The PDT to include
        """
        
        pixel=pdt.pixel
        self.pdts[pixel]=pdt
        self.pixels.append(pixel)

    def decimate(self):
        """
        Method to decimate over the PDTs
        
        Notes
        -----
        This decimation is done in place and cannot be undone, without 
        reloading the table.
        """

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
