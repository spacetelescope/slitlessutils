from .simulated import SimulatedFile
from .observed import ObservedFile
from ....logger import LOGGER

def load_wfss(ftype,*args,**kwargs):
    ''' A factory function to generate different WFSSFile types '''
        
    if ftype=='simulated':
        obj=SimulatedFile(*args,**kwargs)
    elif ftype=='observed':
        obj=ObservedFile(*args,**kwargs)
    else:
        LOGGER.warn(f'File type {ftype} is invalid.')
        obj=None
    return obj
            
            
