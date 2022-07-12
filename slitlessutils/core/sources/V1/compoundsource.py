import numpy as np

from ...logger import LOGGER
from ..photometry import SED,Throughput
from .source import Source

class CompoundSource(Source):
    def __init__(self,img,seg,segid,**kwargs):
        Source.__init__(self,img,seg,segid,**kwargs)

        self.seds={}
        for x,y,w in self:
            xd,yd=self.image_coordinates(x,y,dtype=int)
            self.seds[(self.segid,xd,yd)]=SED()

    def seditems(self):
        for (segid,x,y),sed in self.seds.items():
            yield (sed,x,y)

    def set_sed(self,lamb,flam,flamunc=None,band=None,x=None,y=None,**kwargs):
        key=(self.segid,x,y)
        if key in self.seds:
            self.seds[key].set_sed(lamb,flam,flamunc=flamunc)

            if isinstance(band,Throughput):
                xd,yd=self.image_coordinates((x,y),dtype=int)
                g=np.where((self.x==xd) & (self.y==yd))[0]

                self.seds[key].normalize(band,self.fnu*self.w[g])
        else:
            LOGGER.warning(f'SED Key is not found {key}')

    def make_HDU(self,**kwargs):

        nwav=len(self.extpars)

        # make and fill an output 3d flux cube
        dat=np.full((nwav,self.NAXIS[1],self.NAXIS[0]),np.nan,dtype=float)
        unc=np.full((nwav,self.NAXIS[1],self.NAXIS[0]),np.nan,dtype=float)
        for (segid,x,y),sed in self.seds.items():
            xx=int(x+self.ltv[0])
            yy=int(y+self.ltv[1])

            dat[:,yy,xx]=sed.flam
            unc[:,yy,xx]=sed.flamunc

        # make a header from the astrometry
        hdr=self.to_header()

        # add header keywords to support the 3d cube
        hdr['EXTNAME']=('COMPOUND','extension name')
        hdr['EXTVER']=(self.segid,'extension version number')
        hdr['SRCTYPE']=('COMPOUND','is this a simple or compound source')


        # make an HDU out of the datacube and header
        hdu=fits.ImageHDU(data=dat,header=hdr)
        super().update_HDU(hdu,**kwargs)

        return hdu
