from astropy.io import fits
from astropy.wcs import WCS as astropyWCS,utils
import numpy as np
import os

try:
    from jwst import datamodels
    HAS_JWST=True
except:
    HAS_JWST=False


#pixelarea
#pixelarea_map
#cd



#ad2xy
#xy2ad
#xy2xy



# internal
#bitpix (for mkhdr)
#hassip
#putsip
#getrot
#pixelscale




class WCS(astropyWCS):

    
    def __init__(self,hdr):
        astropyWCS.__init__(self,header=hdr)
        self.naxis1=hdr['NAXIS1']
        self.naxis2=hdr['NAXIS2']
        self.NAXIS=(self.naxis1,self.naxis2)

        # check some things
        if np.isnan(self.wcs.equinox):
            self.wcs.equinox=2000.
            
    
        self.ltv=np.array([hdr.get('LTV1',0.),hdr.get('LTV2',0.)],dtype=float)


    def pixelarea(self,x=None,y=None,unit='arcsec2'):
        # estimate pixel scale using approximate algorithm
        # from https://trs.jpl.nasa.gov/handle/2014/40409
        if x is None:
            x=self.crpix[0]
        if y is None:
            y=self.crpix[1]

        
        l0, phi0 = np.deg2rad(self.xy2ad(x-0.5,y-0.5))
        l1, phi1 = np.deg2rad(self.xy2ad(x-0.5,y+0.5))
        l2, phi2 = np.deg2rad(self.xy2ad(x+0.5,y+0.5))
        l3, phi3 = np.deg2rad(self.xy2ad(x+0.5,y-0.5))
        area=np.abs((l1-l3)*(np.sin(phi0)-np.sin(phi2))+
                    (l0-l2)*(np.sin(phi3)-np.sin(phi1)))/2.

        if unit=='sr':
            area*=42545170296.15221   # = (180/np.pi*3600.)**2
            
            
        return area

    def relative_pixelarea(self,**kwargs):
        # can get this from a reference file
        a0=self.pixelarea(units='arcsec2')
        area=self.pixelarea(**kwargs,units='arcsec2')
        return area/a0

        
    
    def xy2xy(self,xi,yi,wcs):
        a,d=self.xy2ad(xi,yi)
        xo,yo=wcs.ad2xy(a,d)
        return xo,yo


    @staticmethod
    def dtype_to_bitpix(self,dtype):
        ''' return the numpy dtype for a specified BITPIX value '''

        if dtype in (np.float64,float):
            return -64
        elif dtype == np.float32:
            return -32
        elif dtype in (np.int64,int):
            return 64
        elif dtype == np.int32:
            return 32
        elif dtype == np.int16:
            return 16
        else:
            return -64
    
        
class ClassicWCS(WCS):


    def __init__(self,hdr,**kwargs):
        WCS.__init__(self,hdr,**kwargs)


    @classmethod
    def from_fits(cls,filename,ext=1,**kwargs):
        hdr=fits.getheader(filename,ext=ext)
        return cls(hdr,**kwargs)



    def ad2xy(self,a,d):
        if isinstance(a,(list,tuple,np.ndarray)):
            if len(a)==1:
                return self.all_world2pix(a[0],d[0],0)
            else:
                return self.all_world2pix(a,d,0)
        else:
            return self.all_world2pix(a,d,0)
            
    def xy2ad(self,x,y):
        if isinstance(x,(list,tuple,np.ndarray)):
            if len(x)==1:
                return self.all_pix2world(x[0],y[0],0)
            else:
                return self.all_pix2world(x,y,0)
        else:
            return self.all_pix2world(x,y,0)

    #def pixelarea(self,x=None,y=None,unit='arcsec2'):
    #    pass
    #
    #def relative_pixelarea(self,*kwargs):
    #    pass
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

    

    def mkhdr(self,dtype,**kwargs):
        if self.wcs.ctype[0][:8]=='RA---TAN' and \
           self.wcs.ctype[1][:8]=='DEC--TAN':
            print("GET CODE AND AUTHOR")
            hdr=fits.Header()
            hdr['SIMPLE']=(True,'')#f'Created by {__code__} ({__author__})')

            bitpix=self.BITPIX.get(dtype)
            if bitpix is None:
                raise NotImplementedError(f"NumPy type of {dtype} is unknown.")
            hdr['BITPIX']=(bitpix,'Number of bits per data pixel')

            
            now=datetime.datetime.now().strftime("%Y-%m-%d")
            hdr['DATE']=(now,'Creation UTC (CCCC-MM-DD) date of FITS header')
            hdr['NAXIS']=(self.wcs.naxis,'Number of data axes')
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
                self._putSIP(hdr,self.sip.a,'A')
                self._putSIP(hdr,self.sip.b,'B')
                self._putSIP(hdr,self.sip.ap,'AP')
                self._putSIP(hdr,self.sip.bp,'BP')


            # add any additional keywords
            for k,v in kwargs.items():
                hdr[k]=v
       
        else:
            raise NotImplementedError("Invalid CTYPE")

        return hdr



    def _putSIP(self,hdr,sip,name):
        if sip is not None:
            ii,jj=np.where(sip !=0)
            order=getattr(self.sip,f'{name.lower()}_order')
            hdr[f'{name}_ORDER']=order
            for (i,j,v) in zip(ii,jj,sip[ii,jj]):
                hdr[f'{name}_{i}_{j}']=v

    #@staticmethod
    #def _hasSIP(hdr):
    #    keys=list(hdr.keys())        
    #    for func in ('A','B','AP','BP'):
    #        # build a regex search
    #        r = re.compile(f"{func}_[0-9]_[0-9]")
    #        
    #        # has a given order
    #        if f'{func}_ORDER' in keys and any(filter(r.match,keys)):
    #            return True
    #    return False
            



    

class GeneralizedWCS(WCS):
    def __init__(self,fitsfile,ext=1):
        self.fitsfile=fitsfile
        hdr=fits.getheader(self.fitsfile,ext=ext)
        WCS.__init__(self,hdr)
        
        with datamodels.open(self.fitsfile) as dm:
            self._ad2xy=dm.meta.wcs.get_transform('world','detector')
            self._xy2ad=dm.meta.wcs.get_transform('detector','world')


        # pixel scales are in
        #print(self._xy2ad[3].c1_0.value)
        #print(self._xy2ad[4].c0_1.value)

        
        # set the CRPIX
        if self.wcs.crpix[0]==0:
            self.wcs.crpix[0]=self._ad2xy[-3].parameters[0]+1
        if self.wcs.crpix[1]==0:
            self.wcs.crpix[1]=self._ad2xy[-2].parameters[0]+1

        # set the CRVAL
        crval=self.xy2ad(self.wcs.crpix[0],self.wcs.crpix[1])
        if self.wcs.crval[0]==0:
            self.wcs.crval[0]=crval[0]    #+self._ad2xy[1].parameters[0]
        if self.wcs.crval[1]==0:
            self.wcs.crval[1]=crval[1]    #-self._ad2xy[1].parameters[1]

        # set the PC matrix
        if np.array_equal(self.wcs.pc, np.eye(2)):
            pa=-np.deg2rad(self._ad2xy[1].parameters[2])
            cs=np.cos(pa)
            sn=np.sin(pa)
            self.wcs.pc=np.array([[-cs,sn],[sn,cs]])

        # set the CDELT
        if self.wcs.cdelt[0]==1.0:
            self.wcs.cdelt[0]=self._xy2ad[3].c1_0.value/3600.
        if self.wcs.cdelt[1]==1.0:
            self.wcs.cdelt[1]=self._xy2ad[4].c0_1.value/3600.

        # set the CTYPE?

        
    def ad2xy(self,a,d):
        if isinstance(a,(list,tuple,np.ndarray)):
            if len(a)==1:
                x,y,_,_=self._ad2xy(a[0],d[0],0,0)
            else:
                z=np.zeros_like(a,dtype=int)
                x,y,_,_=self._ad2xy(a,d,z,z)
        else:
            x,y,_,_=self._ad2xy(a,d,0,0)
        return x,y

    def xy2ad(self,x,y):
        if isinstance(x,(list,tuple,np.ndarray)):
            if len(x)==1:
                a,d,_,_=self._xy2ad(x[0],y[0],0,0)
            else:
                z=np.zeros_like(x,dtype=int)
                a,d,_,_=self._xy2ad(x,y,z,z)
        else:
            a,d,_,_=self._xy2ad(x,y,0,0)
        return a,d
    


    def mkhdr(self,dtype,**kwargs):
        raise NotImplementedError



def load_wcs(init,**kwargs):
    if isinstance(init,str):
        if os.path.basename(init)=='.fits':
            h=fits.getheader(init,ext=0)
            tel=h.get('TELESCOP','')
            if tel == 'HST':
                obj=ClassicWCS.from_fits(init,**kwargs)
            elif tel=='JWST':
                if h['EXP_TYPE'] in ('NIS_WFSS',):
                    if HAS_JWST:
                        obj=GeneralizedWCS(init,**kwargs)
                    else:
                        print("JWST IS NOT INSTALLED")
                else:
                    print(h['EXP_TYPE'])
                    obj=ClassicWCS.from_fits(init,**kwargs)
            else:
                print("Unknown telescope")
        else:
            print("Unknown file type")
    elif isinstance(init,fits.Header):
        obj=ClassicWCS(init,**kwargs)
    else:
        print("Unknown datatype")

    return obj

        

if __name__=='__main__':

    x=GeneralizedWCS('/Users/rryan/JWST/SMACS0723/MAST_2022-07-13T1120/JWST/jw02736003001_02102_00001_nis/jw02736003001_02102_00001_nis_rate.fits')


