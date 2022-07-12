import numpy as np


def load(h5,key,ptype=None):
    ''' load attributes from the HDF5 file '''
    
    val=h5.attrs.get(key,None)
    if isinstance(val,bytes):
        val=val.decode('UTF-8')
        vlo=val.lower()
        if vlo=='true':
            val=True
        elif vlo=='false':
            val=False
        elif vlo=='none':
            val=None
        else:
            pass
    elif isinstance(val,float):
        if np.isnan(val):
            val=None
        else:
            pass
    else:
        pass
    
    if isinstance(ptype,type):
        val=ptype(val)
            
    return val


def write(h5,key,val):
    ''' write an attribute to the hdf5 object '''

    
    if val is not None:
        if isinstance(val,(bool,str)):
            h5.attrs[key]=np.string_(val)
        else:
            h5.attrs[key]=val


