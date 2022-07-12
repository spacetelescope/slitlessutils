import numpy as np


from .hdf5table import HDF5Table,attributes



class PDT(HDF5Table):
    COLUMNS=('x','y','lam','val')

    def __init__(self,x,y,dims=None,**kwargs):
        HDF5Table.__init__(self,dims=dims,**kwargs)
        self.pixel=(x,y)

        # record the pixel in the attrs
        self.attrs['x']=self.DTYPES['x'](x)
        self.attrs['y']=self.DTYPES['y'](y)
        
    def wavelengths(self,lam=None):

        if 'wav0' in self.attrs and 'dwav' in self.attrs:
            wav0=self.attrs['wav0']
            dwav=self.attrs['dwav']
            
            if lam is None:
                lam=self.get('lam')

            wav=wav0+lam*dwav
            return wav
        else:
            return lam
        
    @property
    def name(self):
        #'({},{})'.format(*self.pixel)
        return str(self.pixel)




    @classmethod
    def load_hdf5(cls,h5,x,y):
        
        obj=cls(x,y)

        # open the dataset
        hd=h5[obj.name]

        # load the attributes
        for attr in hd.attrs:
            obj.attrs[attr]=attributes.load(hd,attr)
        
        # transfer the data into the columns
        data=hd[()]
        for column in obj.COLUMNS:
            obj[column]=data[column]

        return obj
        


