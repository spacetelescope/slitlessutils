import os
from astropy.io import fits
import numpy as np

from .wfss import WFSSFile,WFSSImage
from ....logger import LOGGER

class ObservedImage(WFSSImage):
    TYPE='observed'
    def __init__(self,filename,detconf):
        self.filename=filename
        hdr=fits.getheader(self.filename,
                           extname=detconf.extensions['science'].name,
                           extver=detconf.extensions['science'].ver)

        WFSSImage.__init__(self,hdr,detconf)


    def exten(self,imgtype):
        return (self.extensions[imgtype].name,self.extensions[imgtype].ver)


    def readfits(self,imgtype,countrate=True,header=True):
        ''' basic utility to read fits data from disk '''
        
        exten=(self.extensions[imgtype].name,self.extensions[imgtype].ver)

        img,hdr=fits.getdata(self.filename,exten,header=True)

        if countrate:
            if 'bunit' in hdr:
                bunit=hdr['bunit'].lower()
                if bunit in ('electron','electrons','e','e-',
                             'datanumber','datanumbers','dn'):
                    h0=fits.getheader(self.filename,0)
                    if 'exptime' in h0:
                        img/=h0['exptime']
                    else:
                        LOGGER.warning(f"BUNIT implies count image, but no EXPTIME for {self.filename}")
            else:
                LOGGER.warning(f"BUNIT is missing, cannot check units {self.filename}")
                        

        if header:
            return img,hdr
        else:
            return img

     
    def read_science(self,**kwargs):
        return self.readfits('science',**kwargs)

    def read_uncertainty(self,**kwargs):
        return self.readfits('uncertainty',**kwargs)

    def read_dataquality(self,**kwargs):
        return self.readfits('dataquality',**kwargs)

    def headfits(self,imgtype):
        return fits.getheader(self.filename,self.exten(imgtype))

    
class ObservedFile(WFSSFile):
    TYPE='observed'
    def __init__(self,filename,insconf,**kwargs):

        dataset='_'.join(filename.split('_')[:-1])

        self.dataset=dataset
        self.gzip=dataset[-3:]=='.gz'

        WFSSFile.__init__(self,self.dataset,insconf)

        # read the primary header
        self.phdr=fits.getheader(filename,ext=0)

        for detname,detconf in insconf.items():
            self[detname]=ObservedImage(filename,detconf)
