import numpy as np
import datetime
from astropy.io import fits
from astropy.wcs import WCS as astropyWCS
import astropy.wcs.utils as wcsutils
import re

from ...logger import LOGGER
from ...info import __code__,__author__

def fixscalar(func):
    def wrapper(*args,**kwargs):
        a,b=func(*args,**kwargs)
        if isinstance(a,np.ndarray) and a.ndim==0:
            a=a.item()
            b=b.item()
        return a,b
    return wrapper



class WCS(astropyWCS):
    ''' A class to override the astropy.wcs.WCS to do few more things '''

    BITPIX={np.float64: -64,
            np.float32: -32,
            np.int64: 64,
            np.int32: 32,
            np.int16: 16,
            float: -64,
            int: 64}
   
    
    def __init__(self,hdr):
        # ok, there is this annoying situation where the SIP coefs will
        # be present, but the CTYPE does not indicate that.  Therefore,
        # astropy will issue a big warning, so I'll take control of that
        # here, but checking for SIP coefs and the CTYPE
        #has_sip=self.has_sip(hdr)
        #if has_sip:
        #    for ctype in ('CTYPE1','CTYPE2'):
        #        if hdr[ctype][8:] != '-SIP':
        #            hdr[ctype] += '-SIP'

        # set the astropy WCS class
        astropyWCS.__init__(self,header=hdr)
        
        self.naxis1=hdr['NAXIS1']
        self.naxis2=hdr['NAXIS2']
        
        self.NAXIS=(self.naxis1,self.naxis2)
        #self.shape=(self.naxis2,self.naxis1)
        self.npixel=self.naxis1*self.naxis2
        
        # check some things
        if np.isnan(self.wcs.equinox):
            self.wcs.equinox=2000.
                        
        # for reasons, they do not include LTV* in WCS
        self.ltv=np.array([hdr.get('LTV1',0.),hdr.get('LTV2',0.)],dtype=float)
     
        
    def bitpixToNumpy(self,bitpix):
        ''' return the numpy dtype for a specified BITPIX value '''

        for k,v in self.BITPIX.items():
            if v==bitpix:
                return k
        else:
            return np.float64

        
   
    def pixel_area_map(self,x,y,scale=1.):
        ''' compute the pixel area for the distortion '''
        
        dx=np.array(x)-(self.wcs.crpix[0]-1.)
        dy=np.array(y)-(self.wcs.crpix[1]-1.)
        try:
            n=len(x)
        except:
            n=1

        # compute derivatives for A
        dadx=np.ones_like(dx)
        dady=np.zeros_like(dx)
        if hasattr(self.sip,'a'):
            ii,jj=np.where(self.sip.a != 0)
            for i,j,v in zip(ii,jj,self.sip.a[ii,jj]):
                if i!=0:
                    dadx+=(v*i*dx**(i-1)*dy**j)
                if j!=0:
                    dady+=(v*j*dx**i*dy**(j-1))
                

        # compute derivatives for B         
        dbdx=np.zeros_like(dx)
        dbdy=np.ones_like(dx)
        if hasattr(self.sip,'b'):
            ii,jj=np.where(self.sip.b != 0)
            for i,j,v in zip(ii,jj,self.sip.b[ii,jj]):
                if i!=0:
                    dbdx+=(v*i*dx**(i-1)*dy**j)
                if j!=0:
                    dbdy+=(v*j*dx**i*dy**(j-1))

        jacobian = scale*np.abs(dadx * dbdy - dady * dbdx)
        
        if n==1:
            jacobian=jacobian[0]
        
        return jacobian



    @staticmethod
    def has_sip(hdr):

        keys=list(hdr.keys())        
        for func in ('A','B','AP','BP'):
            # build a regex search
            r = re.compile(f"{func}_[0-9]_[0-9]")

            
            # has a given order
            if f'{func}_ORDER' in keys and any(filter(r.match,keys)):
                return True
        return False
            

            
    def putSIP(self,hdr,sip,name):
        if sip is not None:
            ii,jj=np.where(sip !=0)
            order=getattr(self.sip,f'{name.lower()}_order')
            hdr[f'{name}_ORDER']=order
            for (i,j,v) in zip(ii,jj,sip[ii,jj]):
                hdr[f'{name}_{i}_{j}']=v
        
    def mkhdr(self,dtype,**kwargs):       

        if self.wcs.ctype[0][0:8]=='RA---TAN' and \
           self.wcs.ctype[1][0:8]=='DEC--TAN':

            hdr=fits.Header()
            hdr['SIMPLE']=(True,f'Created by {__code__} ({__author__})')

            bitpix=self.BITPIX.get(dtype)
            if bitpix is None:
                raise NotImplementedError(f"NumPy type of {dtype} is unknown.")
            hdr['BITPIX']=(bitpix,'Number of bits per data pixel')

            
            now=datetime.datetime.now().strftime("%Y-%m-%d")
            hdr['DATE']=(now,'Creation UTC (CCCC-MM-DD) date of FITS header')
            hdr['NAXIS']=(2,'Number of data axes')
            hdr['NAXIS1']=(self.shape[1],'Number of pixels in x-axis')
            hdr['NAXIS2']=(self.shape[0],'Number of pixels in y-axis')
            
            # put in the astrometry
            hdr['CRPIX1']=(self.wcs.crpix[0],'x-coordinate of reference pixel')
            hdr['CRPIX2']=(self.wcs.crpix[1],'y-coordinate of reference pixel')
            hdr['CRVAL1']=(self.wcs.crval[0],'first axis value at reference pixel')
            hdr['CRVAL2']=(self.wcs.crval[1],'second axis value at reference pixel')

            # astropy complains about this, ugh.
            #hdr['CDELT1']=(self.wcs.cdelt[0],' ')
            #hdr['CDELT2']=(self.wcs.cdelt[1],' ')

            hdr['CTYPE1']=(self.wcs.ctype[0],'the coordinate type for the first axis') 
            hdr['CTYPE2']=(self.wcs.ctype[1],'the coordinate type for the second axis')

            cd=self.cd(unit='deg')
            hdr['CD1_1']=(cd[0,0],'partial of first axis coordinate w.r.t. x') 
            hdr['CD2_1']=(cd[1,0],'partial of second axis coordinate w.r.t. x')
            hdr['CD1_2']=(cd[0,1],'partial of first axis coordinate w.r.t. y')
            hdr['CD2_2']=(cd[1,1],'partial of second axis coordinate w.r.t. y')

            #if hasattr(self,'ltv'):
            hdr['LTV1']=(self.ltv[0],'x-axis pixel reference')
            hdr['LTV2']=(self.ltv[1],'y-axis pixel reference')

            if np.isnan(self.wcs.equinox):
                hdr['EQUINOX']=(2000.,'Equinox of Ref. Coord.')
            else:
                hdr['EQUINOX']=(self.wcs.equinox,'Equinox of Ref. Coord.')
            hdr['LONGPOLE']=(self.wcs.lonpole,' ')

            
            if self.wcs.ctype[0][9:]==self.wcs.ctype[1][9:]=='SIP':
                self.putSIP(hdr,self.sip.a,'A')
                self.putSIP(hdr,self.sip.b,'B')
                self.putSIP(hdr,self.sip.ap,'AP')
                self.putSIP(hdr,self.sip.bp,'BP')


            # add any additional keywords
            for k,v in kwargs.items():
                hdr[k]=v
       


        else:
            raise NotImplementedError("Invalid CTYPE")

        return hdr


    def cd(self,unit='deg'):

        try:
            cd=self.wcs.cd
        except:
            #if ~self.wcs.has_cd():
            cd=self.pixel_scale_matrix
            cd[0,:]*=self.wcs.cdelt[0]
            cd[1,:]*=self.wcs.cdelt[1]            

 

        if unit=='rad':
            cd=np.deg2rad(cd)
        elif unit=='arcsec':
            cd*=3600.
        elif unit=='arcmin':
            cd*=60.
        elif unit=='deg':
            pass            
        else:
            raise NotImplementedError(f'invalud unit: {unit}')

        return cd
        

    
    def getrot(self):

        cd=self.cd(unit='rad')

        # astropy throws this as a warning
        #if self.wcs.cdelt[0] != 1:
        #    cd[0,:]=cd[0,:]*self.wcs.cdelt[0]
        #if self.wcs.cdelt[1] != 1:
        #    cd[1,:]=cd[1,:]*self.wcs.cdelt[1]
        

            
        det=cd[0,0]*cd[1,1]-cd[0,1]*cd[1,0]
        if det < 0:
            sgn=-1
        else:
            sgn=+1

        if (cd[1,0] == 0) and (cd[0,1]==0):
            rot=0.
            pix=np.array([cd[0,0],cd[1,1]])
        else:
            rot1=np.arctan2(sgn*cd[0,1],sgn*cd[0,0])
            rot2=np.arctan2(-cd[1,0],cd[1,1])

            
            if rot1 != rot2:
                if np.rad2deg(np.abs(rot1-rot2))<2:
                    rot=(rot1+rot2)/2.
                elif np.rad2deg(np.abs(rot1-rot2-2*np.pi))<2:
                    rot=(rot1+rot2-2*np.pi)/2
                elif np.rad2deg(np.abs(rot1-rot2+2*np.pi))<2:
                    rot=(rot1+rot2+2*np.pi)/2
                else:
                    LOGGER.warn("x/y axis rotations differ by >=2 deg")
                    rot=rot1
            else:
                rot=rot1
                    
                    
            
            pix=np.array([sgn*np.sqrt(cd[0,0]*cd[0,0]+cd[0,1]*cd[0,1]),\
                          np.sqrt(cd[1,1]*cd[1,1]+cd[1,0]*cd[1,0])])
            
            rot=np.rad2deg(rot)
            
            if self.wcs.lonpole != 180:
                rot=rot+(180-self.wcs.lonpole)

        pix=np.rad2deg(pix)
        return rot,pix

    @property
    def pixelscale(self):
        ''' Return pixel scale in arcsec '''
        
        #rot,cdelt=self.getrot()
        #scales=np.abs(cdelt)
        scales=wcsutils.proj_plane_pixel_scales(self)        
        return scales*3600.

        
        
    @property
    def pixelarea(self):
        ''' Return pixel area in arcsec2 '''
        #rot,cdelt=self.getrot()
        #area=-(cdelt[0]*3600)*(cdelt[1]*3600)
        pscl=self.pixelscale
        area=pscl[0]*pscl[1]        
        return area


    # convenience functions because that zero is a nuissance
    @fixscalar
    def xy2ad(self,x,y):
        return self.all_pix2world(x,y,0)

    @fixscalar
    def ad2xy(self,a,d):
        return self.all_world2pix(a,d,0)

    @fixscalar
    def xy2xy(self,x,y,obj):
        assert isinstance(obj,astropyWCS)
        a,d=self.all_pix2world(x,y,0)
        return obj.all_world2pix(a,d,0)
        
#looking at pickling?    
#    def __reduce__(self):
#        """
#        Support pickling of WCS objects.  This is done by serializing
#        to an in-memory FITS file and dumping that as a string.
#        """
#
#        hdulist = self.to_fits(relax=True)
#        
#        
#        buffer = io.BytesIO()
#        hdulist.writeto(buffer)
#
#        dct = self.__dict__.copy()
#        dct['_alt_wcskey'] = self.wcs.alt
#        
#        out=(__unpickle__,(self.__class__, dct, buffer.getvalue(),))
#        
#        return out
#             
#
#
#    
#def __unpickle__(cls, dct, fits_data):
#    """
#    Unpickles a WCS object from a serialized FITS string.
#    """
#
#    print('unpickle')
#    self = cls.__new__(cls)
#    
#    buffer = io.BytesIO(fits_data)
#    hdulist = fits.open(buffer)
#    
#    naxis = dct.pop('naxis', None)
#    if naxis:
#        hdulist[0].header['naxis'] = naxis
#        naxes = dct.pop('_naxis', [])
#        for k, na in enumerate(naxes):
#            hdulist[0].header[f'naxis{k + 1:d}'] = na
#            
#    kwargs = dct.pop('_init_kwargs', {})
#    self.__dict__.update(dct)
#    
#    wcskey = dct.pop('_alt_wcskey', ' ')
#    WCS.__init__(self, hdulist[0].header, hdulist, key=wcskey, **kwargs)
#    self.pixel_bounds = dct.get('_pixel_bounds', None)
#    
#    return self

        

        
if __name__=='__main__':
    
    FILE='/Users/rryan/icoi3immq_flt.fits'


    import multiprocessing as mp
    import pickle

    
    with fits.open(FILE) as hdul:
        hdr=hdul[1].header
    wcs=WCS(hdr)
        
    xys=[(i,i+1) for i in range(4)]
    it=((wcs,xy) for xy in xys)

    a,d=wcs.xy2ad(500,500)
    print(a,d)

    
    #ads=[wcs.xy2ad(x,y) for wcs,(x,y) in it]
    with open('t.pickle','ab') as f:
        pickle.dump(wcs,f)
    with open('t.pickle','rb') as f:
        ww=pickle.load(f)

    a,d=ww.xy2ad(500,500)
    a,d=ww.xy2ad(500,500)
    print(a,d)
    

    ljkj
    
    with mp.Pool(processes=2) as p:
        imap=p.imap(test,it)
        ads=list(imap)

