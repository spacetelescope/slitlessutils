import numpy as np
from astropy.io import fits

from ...utilities import headers


def check_in_field(func):
    ''' decorator to check if (x,y) is within image bounds '''
    def wrapper(self,x,y,l):

        xx=np.rint(x).astype(int)
        yy=np.rint(y).astype(int)

        logic=(xx > 0) & (xx < self.shape[1]) & (yy > 0) & (yy < self.shape[0])
        f=np.where(logic,func(self,xx,yy,l),0.)
        return f.astype(self.DTYPE)
    return wrapper



class FlatField:
    ''' base class for all flat field types '''

    DTYPE = np.float32

    def __init__(self,filename):
        self.filename=filename


    def update_header(self,hdr):
        hdr['FFTYPE']=(self.FTYPE,'flat field image type')
        hdr['FFNAME']=(self.filename,'flat field image name')
        headers.add_stanza(hdr,'Flat Field Settings',before='FFTYPE')





class UnityFlatField(FlatField):

    ''' Unity flat field '''

    FTYPE = 'unity'

    def __init__(self):
        FlatField.__init__(self,'')


    def __call__(self,x,y,l):
        return np.ones_like(x,dtype=self.DTYPE)

    def update_header(self,hdr):
        super().update_header(hdr)


class ImageFlatField(FlatField):
    ''' Image/Gray flat field (ie. no wavelength dependence) '''


    FTYPE='gray'


    def __init__(self,filename):
        FlatField.__init__(self,filename)

    @classmethod
    def from_fits(cls,filename):
        obj=cls(filename)


        try:
            data=fits.getdata(obj.filename,ext=1)
            obj.shape=obj.data.shape
        except:
            LOGGER.warning("Image flat is invalid, using unity")
            obj=UnityFlatField()


        return obj

    @check_in_field
    def __call__(self,x,y,l):
        return self.data[y,x]

    def update_header(self,hdr):
        super().update_header(hdr)


class PolynomialFlatField(FlatField):
    ''' Flatfield with polynomial-wavelength dependence '''

    FTYPE='polynomial'
    def __init__(self,filename):
        FlatField.__init__(self,filename)

    @classmethod
    def from_fits(cls,filename):
        obj=cls(filename)

        with fits.open(obj.filename,mode='readonly') as hdul:
            obj.data=[hdu.data for hdu in hdul if hdu.data is not None]
        obj.order=len(obj.data)-1


        if obj.order ==-1:
            LOGGER.warning("Polynomial flat is invalid, using unity")
            obj=UnityFlatField()
        elif obj.order==0:
            obj=ImageFlatField.from_fits(filename)
        else:
            obj.wmin=hdul[0].header['WMIN']
            obj.wmax=hdul[0].header['WMAX']
            obj.shape=obj.data[0].shape

        return obj

    @check_in_field
    def __call__(self,x,y,l):
        ll=(l-self.wmin)/(self.wmax-self.wmin)
        return sum(f[y,x]*ll**i for i,f in enumerate(self.data))

    def update_header(self,hdr):
        super().update_header(hdr)



def load_flatfield(*args,unity=False):
    ''' A factory function to load different flat field types '''

    n=len(args)
    if n==0 or unity:
        obj=UnityFlatField()
    elif n==1:
        if isinstance(args[0],str):
            obj=PolynomialFlatField.from_fits(args[0])
        else:
            raise NotImplementedError(f"data type of args: {type(args)}")
    else:
        raise NotImplementedError(f"Unknown number of arguments")
    return obj
