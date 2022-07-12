import numpy as np
from astropy.io import fits


from .wcs import WCS
from ..photometry import Throughput
from ...logger import LOGGER
class AstroImage(WCS):


    def __init__(self,image,header,zeropoint=None,throughput=None):
        self.image=image
        self.header=header
        WCS.__init__(self,self.header)


        # load somethings
        self.zeropoint=self._get_zeropoint() if zeropoint is None else zeropoint
        self.throughput=self._get_throughput() if throughput is None else throughput



    def _get_zeropoint(self,default=25.0):
        ''' set the zeropoint from the header '''
        zeropoint=None

        bunit=self.header.get('BUNIT','').lower()
        if bunit=='mjy/sr':
            self.zeropoint=-2.5*np.log10((1e6/3631.)*(np.pi/180)**2*(self.pixelscale/3600)**2)
        elif bunit in ('electrons/s','electron/s','e/s','e-/s'):
            for k in ('ABZERO','ZEROPT','MAGZERO','ZERO','ZPT'):
                if k in self.header:
                    zeropoint=self.header[k]
                else:
                    # if here, then no zeropoint is found
                    pass
        else:
            # if here, then unknown image units, which is ok (like segmap)
            pass


        if zeropoint is None:
            LOGGER.warning(f"No zeropoint was loaded, using {default}")
            zeropoint=default
        return zeropoint


    def _get_throughput(self):
        ''' load a throughput curve '''

        if 'FILTFILE' in self.header:
            thru=Throughput.from_file(self.header['FILTFILE'])
        elif all(k in self.header for k in ('TELESCOP','INSTRUME','FILTER')):
            thru=Throughput.from_keys(self.header['TELESCOP'],
                                      self.header['INSTRUME'],
                                      self.header['FILTER'])
        else:
            LOGGER.warning("Unable to load any filter curve")
            thru=None

        return thru


    def __eq__(self,a):
        return self.image ==a

    def __ne__(self,a):
        return self.image !=a

    def __ge__(self,a):
        return self.image >=a

    def __gt__(self,a):
        return self.image >a

    def __lt__(self,a):
        return self.image <a

    def __le__(self,a):
        return self.image <=a

    def __sub__(self,a):
        self.image-= a
        return self

    def __isub__(self,a):
        return self.__sub__(a)

    def __rsub__(self,a):
        self.image = a-self.image
        return self

    def __add__(self,a):
        self.image+=a
        return self

    def __radd__(self,a):
        return self.__add__(a)

    def __iadd__(self,a):
        return self.__add__(a)

    def __mul__(self,a):
        self.image*=a
        return self

    def __rmul__(self,a):
        return self.__mul__(a)

    def __imul__(self,a):
        return self.__mul__(a)

    def __contains__(self,k):
        return k in self.header

    def get(self,*args):
        return self.header.get(*args)

    @property
    def shape(self):
        return self.image.shape

    def astype(self,dtype):
        hdr=self.header.copy()
        hdr['BITPIX']=self.BITPIX[dtype]
        return type(self)(self.image.astype(dtype),hdr)


    def extract(self,x0,x1,y0,y1,**kwargs):
        assert (x1>x0 and y1>y0),'Box invalid size'

        # get a bounding box
        ny,nx=self.image.shape
        xx0,xx1=max(x0,0),min(x1+1,nx-1)
        yy0,yy1=max(y0,0),min(y1+1,ny-1)

        # cut the image out
        image=self.image[yy0:yy1,xx0:xx1]

        # copy the header over
        header=self.header.copy()
        header['NAXIS1']=image.shape[1]
        header['NAXIS2']=image.shape[0]
        header['CRPIX1']-=xx0
        header['CRPIX2']-=yy0
        header['LTV1']=(self.get('LTV1',0.)-xx0,'difference between phys/img coord')
        header['LTV2']=(self.get('LTV2',0.)-yy0,'difference between phys/img coord')

        # add new keywords
        header['XMIN']=(xx0,'lower x-bound')
        header['XMAX']=(xx1,'upper x-bound')
        header['YMIN']=(yy0,'lower y-bound')
        header['YMAX']=(yy1,'upper y-bound')

        # update the history
        history=f'Extracted from region (x,y)=[{xx0}:{xx1},{yy0}:{yy1}]'
        header.add_history(history)

        # package the output as the same type as the parent
        out=type(self)(image,header,zeropoint=self.zeropoint,throughput=self.throughput)

        return out

    def writefits(self,filename,**kwargs):
        fits.writeto(filename,self.image,self.header,**kwargs)

    @classmethod
    def from_HDU(cls,hdu,**kwargs):
        return cls(hdu.data,hdu.header,**kwargs)


    @classmethod
    def from_fits(cls,filename,exten,**kwargs):
        with fits.open(filename,mode='readonly') as hdul:
            obj=cls.from_HDU(hdul[exten],**kwargs)

        return obj
