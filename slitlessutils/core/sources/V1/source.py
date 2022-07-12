from astropy.io import fits
from skimage.segmentation import expand_labels
import numpy as np
import os

from ..astrometry import WCS
from ..utilities import headers
from ...logger import LOGGER

class Source(WCS):
    def __init__(self,img,seg,segid,mode='pixels',backsize=5,zeropoint=0.0):
        self.segid=segid
        self.mode=mode
        self.backsize=backsize
        self.zeropoint=zeropoint

        LOGGER.debug("must set EXTPARS")


        # compute weights
        self.set_weights(img,seg)

        # compute the center-of-mass in the subimage
        xyc=(np.average(self.x,weights=self.w),
             np.average(self.y,weights=self.w))

        # translate to pixel position
        self.xyc=self.image_coordinates(*xyc)
        self.adc=self.xy2ad(*xyc)
        self.flux=np.sum(img[self.y,self.x])
        self.fnu=self.flux*10.**(-0.4*(self.zeropoint+48.6))
        self.mag=-2.5*np.log10(self.flux)+self.zeropoint

    def __str__(self):
        return f'{self.__class__.__name__} source: {self.segid}'

    def __len__(self):
        return len(self.x)

    def __iter__(self):
        yield from zip(self.x,self.y,self.w)

    def __bool__(self):
        return len(self)>0

    @property
    def name(self):
        return str(self.segid)

    @property
    def xd(self):
        return (self.x-self.ltv[0]).astype(int)

    @property
    def yd(self):
        return (self.y-self.ltv[1]).astype(int)

    def unknowns(self):
        return {key:len(sed) for key,sed in self.seds.items()}

    
    def set_weights(self,img,seg):

        y,x=np.where((img>0) & (seg==self.segid))
        npix=len(x)

        if npix>0:

            # compute a local background
            if self.backsize>0:

                # get a bounding box
                x0=max(np.amin(x)-self.backsize,0)
                y0=max(np.amin(y)-self.backsize,0)
                x1=min(np.amax(x)+self.backsize,img.shape[1]-1)
                y1=min(np.amax(y)+self.backsize,img.shape[0]-1)

                # cut out regions
                subimg=img[y0:y1,x0:x1]
                subseg=seg[y0:y1,x0:x1]

                # grow the segmapsize the segmap
                subseg2=expand_labels(subseg,distance=self.backsize)
                g=np.where((subseg2==segid) & (subseg==0))

                # compute average in background region
                ave,med,sig=sigma_clipped_stats(subimg[g],
                                                sigma_lower=nsiglo,
                                                sigma_upper=nsighi)
                self.background=ave
            else:
                self.background=0.0


            # subtract the background
            img-=self.background

            # process each option
            if self.mode=='pixels':
                self.x=x
                self.y=y
                self.w=np.abs(img[self.y,self.x])

            






            else:
                self.mode='pixels':
                self.x=x
                self.y=y
                self.w=np.abs(img[self.y,self.x])


            self.norm=np.sum(img[self.y,self.x])
            self.w/=self.norm

        else:
            self.x=[]
            self.y=[]
            self.w=[]


    def image_coordinates(self,x,y,dtype=None):
        xx=x-self.ltv[0]
        yy=y-self.ltv[1]
        if dtype is not None:
            xx=dtype(xx)
            yy=dtype(yy)
        return xx,yy

    def instrumental_flux(self,img):
        xd,yd=self.image_coordinates(self.x,self.y,dtype=np.uint32)
        return np.sum(img[yd,xd])

    def write_seds(self,path='.'):
        for sedkey,sed in self.seds.items():
            if isisntance(sedkey,(tuple,list)):
                base='_'.join(str(s) for s in sedkey)
            else:
                base=sedkey
            filename=os.path.join(path,f'{base}.sed')
            sed.write_file(filename)

    def update_HDU(self,hdu):
        hdu.header['LTV1']=self.ltv[0]
        hdu.header['LTV2']=self.ltv[1]

        hdu.header['SEGID']=(self.segid,'Segmentation ID')
        hdu.header['GRPID']=(group,'Group ID')
        hdu.header['RA']=(self.adc[0],'Right Ascension (deg)')
        hdu.header['DEC']=(self.adc[1],'Declination (deg)')
        hdu.header['X']=(self.xyc[0],'X barycenter')
        hdu.header['Y']=(self.xyc[1],'Y barycenter')
        hdu.header['MAG']=(self.mag,'magnitude')
        hdu.header['FLUX']=(self.flux,'instrumental flux in direct image')
        hdu.header['NPIX']=(len(self),'number of extracted pixels')
        hdu.header['AREA']=(self.area,'source area (arcsec2)')
        headers.add_stanza(hdu.header,'Source Properties',before='SEGID')

        for k,v in kwargs.items():
            hdu.header[k]=v

        return hdu


if __name__=='__main__':
    img=fits.getdata('/Users/rryan/Python/Russell-Ryan/slitlessutils/test/gal.fits')
    seg=fits.getdata('/Users/rryan/Python/Russell-Ryan/slitlessutils/test/t_seg.fits')


    src=CompoundSource(img,seg,1)
