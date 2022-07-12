
from ..photometry import SED,Throughput
from .source import Source

class SimpleSource(Source):
    def __init__(self,img,seg,segid,**kwargs):
        Source.__init__(self,img,seg,segid,**kwargs)

        self.seds={self.segid:SED()}


    def seditems(self):
        sed=self.seds[self.segid]
        xd,yd=self.image_coordinates(self.x,self.y,dtype=np.int64)

        for x,y in zip(xd,yd):
            yield (sed,x,y)

    def set_sed(self,lamb,flam,flamunc=None,band=None,**kwargs):
        self.seds[self.segid].set_sed(lamb,flam,flamunc=flamunc)

        # if there is a band passed, then
        if isinstance(band,Throughput):
            self.seds[self.segid].normalize(band,self.fnu)

    def make_HDU(self,**kwargs):
        hdu=self.seds[self.segid].make_HDU()

        hdu.header['EXTNAME']=('SIMPLE','extension name')
        hdu.header['EXTVER']=(self.segid,'extension version number')
        hdu.header['SRCTYPE']=('SIMPLE','is this a simple or composite source')

        super().update_HDU(hdu,**kwargs)

        return hdu
