import numpy as np
from scipy.constants import c
from astropy.io import fits
import os


from .avefnu import avefnu
from ...logger import LOGGER


class SED:
    DTYPE = [('lamb',np.float32),
             ('flam',np.float32),
             ('flamunc',np.float32)]

    def __init__(self,*args,**kwargs):#lamb,flam,flamunc=None):
        #self.set_sed(lamb,flam,flamunc=flamunc)
        if len(args)==2:
            self.set_sed(*args,**kwargs)

    def __bool__(self):
        return hasattr(self,'_data') and len(self._data)>0

    def set_sed(self,lamb,flam,flamunc=None):
        ''' method to set the SED '''
        n=len(lamb)
        if n!= len(flam):
            LOGGER.error('flam and lamb have different lengths')
            return

        if flamunc and (n !=len(flamunc)):
            LOGGER.warning("flamunc is ill-shaped.  using None")
            flamunc=None

        # store results in a np array
        self._data=np.empty(n,dtype=self.DTYPE)
        self._data['lamb']=lamb
        self._data['flam']=flam
        self._data['flamunc']=flamunc
        self.normalizer=1.0


        # get some ranges
        self.wmin=np.amin(lamb)
        self.wmax=np.amax(lamb)

    # a bunch of methods to get the data back
    @property
    def lamb(self):
        return self._data['lamb']

    @property
    def flam(self):
        return self._data['flam']

    @property
    def flamunc(self):
        return self._data['flamunc']

    @property
    def fnu(self):
        return self.flam*(self.lamb/c)*(self.lamb/1e10)

    def __getitem__(self,k):
        return self._data[k]

    # basically implement a shortcut to interpolation
    def __call__(self,wave,fnu=False,**kwargs):
        g=np.where(np.isfinite(self.flam))[0]
        flux=np.interp(wave,self.lamb[g],self.flam[g],**kwargs,
                       left=self.flam[g[0]],right=self.flam[g[-1]])
        if fnu:
            flux*=((wave/c)*(wave/1e10))

        return flux

    def __mul__(self,a):
        self._data['flam']*=a
        self._data['flamunc']*=a
        self.normalizer*=a
        return self

    def __rmul__(self,a):
        return self.__mul__(a)

    def __imul__(self,a):
        return self.__mul__(a)

    def __len__(self):
        try:
            n=len(self._data)
        except:
            n=0
        return n

    def redshift(self,z):
        self._data['lamb']*=(1+z)

    @classmethod
    def from_HDU(cls,hdu):
        # parse data from the HDU
        if 'uncertainty' in hdu.data.names:
            func=hdu.data['uncertainty']
        else:
            func=None
        lamb,flux=hdu.data['wavelength'],hdu.data['flux']

        # adjust units of wavelengths
        units=hdu.header.get('TUNIT1','').lower()
        if units in ('angstrom','angstroms','a',''):
            pass
        elif units in ('micron','microns','um'):
            lamb*=1e4
        else:
            LOGGER.warning(f'Unknown wavelength units {units}. Assuming A')


        # adjust units of flux density
        units=hdu.header.get("TUNIT2",'').lower()
        if units in ('fnu',):
            const=(c/lamb)*(1e10/lamb)
            flam = flux*const
            flamunc = func*const
        elif units in ('flam',):
            flam=np.copy(flux)
            flamunc=np.copy(func)
        else:
            flam=np.copy(flux)
            flamunc=np.copy(func)
            LOGGER.warning(f'Unknown flux units {units}.  Assuming Flam')


        obj=cls(lamb,flam,flamunc=flamunc)

        # check some things in the header
        if 'REDSHIFT' in hdu.header:
            obj.redshift(hdu.header['REDSHIFT'])


        return obj



    @classmethod
    def from_file(cls,filename):
        ''' load an SED from a file (fits, csv, dat, ...) '''

        if os.path.exists(filename):
            ext=os.path.splitext(filename)[1][1:].lower()


            if ext in ('dat','txt','ascii','sed'):
                lamb,flam=np.loadtxt(filename,usecols=(0,1),unpack=True)
                obj=cls(lamb,flam)


            elif ext in ('fits',):
                with fits.open(filename,mode='readonly') as hdul:
                    obj=cls.from_HDU(hdul[1])

                return obj
            elif ext in ('csv',):
                raise NotImplementedError

            else:
                LOGGER.warning(f"File type {ext} is not supported")
                return None
        else:
            LOGGER.warning(f"SED File not found: {filename}")
            return None

        obj.filename=filename      # record the filename

        return obj

    def normalize(self,band,fnu):
        ''' force the flux normalization through a bandpass '''

        factor=fnu/avefnu(self,band)
        self *= factor

        # import matplotlib.pyplot as plt
        # fig,ax=plt.subplots(1,1)
        # ax.plot(self.lamb,self.flam)
        # ax.set_xlim(5000,10000)
        # plt.tight_layout()
        # plt.show()



    def write_file(self,filename,**kwargs):
        ''' write a file (fits, csv, dat, ...) '''


        ext=os.path.splitext(filename)[1][1:].lower()
        if ext in ('csv',):
            delimiter=','
            np.savetxt(filename,self._data,delimiter=delimiter,header=delimiter.join(self._data.dtype.names))
        elif ext in ('txt','dat',''):
            delimiter=' '
            np.savetxt(filename,self._data,delimiter=delimiter,header=delimiter.join(self._data.dtype.names))

        elif ext in ('fits','fit'):
            hdu=self.make_hdu(**kwargs)
            hdu.writeto(filename,overwrite=True)
        else:
            delimiter=' '
            np.savetxt(filename,self._data,delimiter=delimiter,header=delimiter.join(self._data.dtype.names))


    def make_HDU(self,**kwargs):
        ''' package the spectrum into a fits.BinTableHDU '''

        hdu=fits.BinTableHDU(self._data)


        for k,v in kwargs.items():
            hdu.header[k]=v

        return hdu
