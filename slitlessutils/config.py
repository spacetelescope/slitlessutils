import os
import glob

from .logger import LOGGER



# file suffixes.  Probably shouldn't ever change these, but here they are:
SUFFIXES={'1d spectra':'x1d',
          '2d spectra':'x2d',
          'group':'grp',
          'matrix':'mat',
          'wfssimage':'flt'}


class Config:
    ''' Establishes the global configuration '''

    # make this a singleton
    _instance = None
    def __new__(cls,*args,**kwargs):
        if cls._instance is None:
            obj = super().__new__(cls)
            obj.fluxscale = 1e-17               # output numerical flux scale
            obj.fluxunits = 'erg/(s*cm**2*AA)'  # output flux units
            obj.confpath = '/Users/rryan/slitlessutils_config/'
            obj.compression = 'gzip'            # compression type for H5PY
            obj.compression_opts = 9            # compression level for H5PY


                          

            # save the object
            cls._instance=obj

        # set user supplied values
        for k,v in kwargs.items():
            if hasattr(self,k):
                setattr(cls._instance,k,v)
        return cls._instance


    def get_reffile(self,conffile,path='',**kwargs):
        test=os.path.join(self.confpath,path,conffile)
        files=glob.glob(test,**kwargs)

        n=len(files)
        if n==0:
            LOGGER.warn("no config file found")
            return None
        elif n==1:
            return files[0]
        else:
            print("MULTIPLE")



    @property
    def h5pyargs(self):
        return {'compression':self.compression,\
                'compression_opts':self.compression_opts}

    def update_header(self,hdr):
        hdr.set('FLUXSCL',value=self.fluxscale,comment='numeric scale of flux')
        hdr.set('FLUXUNIT',value=self.fluxunits,comment='physical units of flux')
        hdr.set('CONFPATH',value=self.config_path,comment='path to config files')
        hdr.set('COMPTYPE',value=self.compression,comment='type of compression to hdf5')
        hdr.set('COMPOPTS',value=self.compression_opts,comment='options to hdf5 compression')
        hdr.set('',value='',before='FSCALE')
        hdr.set('',value=f'      / Config Settings',before='FSCALE')
        hdr.set('',value='',before='FSCALE')

    def __str__(self):
        pass

if __name__ == '__main__':
    x=Config()
    print(x.fluxscale)
    x.fluxscale=1e-16
    y=Config(fluxscale=2)
    print(repr(y))
    print(x.fluxscale)
    x.compargs['compression']=None
    print(x.compargs)
