import h5py
import os
import numpy as np

from ...logger import LOGGER
from . import attributes
from .hdf5columns import HDF5Columns


#r	Readonly, file must exist (default)
#r+	Read/write, file must exist
#w	Create file, truncate if exists
#w- or x	Create file, fail if exists
#a	Read/write if exists, create otherwise

class HDF5File(HDF5Columns):
    EXT='h5'

    def __init__(self,wfssfile,path='tables',mode='a',remake=True):
        self.dataset=wfssfile.dataset
        self.telescope=wfssfile.telescope
        self.instrument=wfssfile.instrument
        self.grating=wfssfile.grating
        self.pupil='' if wfssfile.pupil is None else wfssfile.pupil
        self.wfsstype=wfssfile.TYPE
        self.remake=remake

                
        self.attrs={}
        

        self.path=path
        self.mode=mode


    @property
    def filename(self):
        base=os.path.join(self.path,self.dataset)
        return f'{base}_{self.TTYPE}.{self.EXT}'
        
    def __enter__(self):
        if self.mode=='w' and os.path.exists(self.filename) and not self.remake:
            return self.filename

        
        self.h5file=h5py.File(self.filename,self.mode)

        if self.mode != 'r':
            attributes.write(self.h5file,'telescope',self.telescope)
            attributes.write(self.h5file,'instrument',self.instrument)
            attributes.write(self.h5file,'grating',self.grating)
            attributes.write(self.h5file,'pupil',self.pupil)
            attributes.write(self.h5file,'wfsstype',self.wfsstype)

        else:
            self.telescope=attributes.load(self.h5file,'telescope')
            self.instrument=attributes.load(self.h5file,'instrument')
            self.grating=attributes.load(self.h5file,'grating')
            self.pupil=attributes.load(self.h5file,'pupil')
            self.wfsstype=attributes.load(self.h5file,'wfsstype')
            
            
        return self        


    def __str__(self):
        return 'HDF5File'
    
    def __exit__(self,et,ev,etb):
        self.close()

    def __del__(self):
        self.close()


    def close(self):
        if hasattr(self,'h5file'):
            if self.h5file:
                if self.mode != 'r':
                    for k,v in self.attrs.items():
                        attributes.write(self.h5file,k,v)
                self.h5file.close()
        else:
            LOGGER.warn('no HDF5file is open to close')
                
        
    def add_detector(self,detconf):
        if hasattr(self,'h5file'):
            self.h5detector=self.h5file.require_group(detconf.name)

            if self.mode != 'r':
                attributes.write(self.h5detector,'naxis1',
                                 self.DTYPES['x'](detconf.naxis[0]))
                attributes.write(self.h5detector,'naxis2',
                                 self.DTYPES['y'](detconf.naxis[1]))
                
        else:
            LOGGER.warn("no HDF5File is loaded")



    def add_order(self,ordconf):
        if hasattr(self,'h5detector'):

            self.h5order=self.h5detector.require_group(ordconf.order)

            if self.mode != 'r':
                # put some stuff here
                attributes.write(self.h5order,'conffile',ordconf.conffile)
                attributes.write(self.h5order,'order',ordconf.order)
    
        else:
            LOGGER.warn("no HDF5Detector is loaded in HDF5File")


    def load_detector(self,detname):
        try:
            self.h5detector=self.h5file[detname]
        except:
            LOGGER.warn(f'HDF5File does not have detector {detname}')


    def load_order(self,ordname):
        try:
            self.h5order=self.h5detector[ordname]
        except:
            LOGGER.warn(f'HDF5Detector does not have order {ordname}')

            
