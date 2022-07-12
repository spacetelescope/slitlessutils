import os

from ..utilities import Pool,headers
from ...logger import LOGGER

class Module:
    def __init__(self,func,path='tables',ncpu=None,postfunc=None,**kwargs):
        self.path=path
        self.ncpu=ncpu
        self.func=func
        self.postfunc=postfunc

    def __str__(self):
        return '\n'.join([f'  path: {self.path}',f'  ncpu: {self.ncpu}'])



    def __call__(self,wfssdata,sources,**kwargs):
        if len(wfssdata)==0:
            LOGGER.warning('there are no WFSS images to process')
            return

        if len(sources)==0:
            LOGGER.warning("There are no sources to process")
            return


        if hasattr(self,'DESCRIPTION'):
            desc=self.DESCRIPTION
        else:
            desc=''

        # build a custom pool object:
        pool=Pool(self.func,ncpu=self.ncpu,desc=desc)


        # process each wfssdatum
        results=pool(wfssdata,sources,total=wfssdata.nfiles,**kwargs)

        # if necessary postprocess
        if callable(self.postfunc):
            results=self.postfunc(results,wfssdata,sources,**kwargs)



        return results

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self,path):
        try:
            if not os.path.isdir(path):
                os.mkdir(path)

            self._path=path
        except:
            LOGGER.warn(f'Cannot make directory: {path}')



    def update_header(self,hdr):
        if not self.ncpu:
            ncpu=1
        else:
            ncpu=self.ncpu

        hdr['NCPU']=(ncpu,'Number of CPUs used')
        hdr['HDF5PATH']=(self.path,'path to HDF5 tables')
        headers.add_stanza(hdr,'Module Settings',before='NCPU')
