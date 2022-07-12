import numpy as np
import os
from skimage.segmentation import expand_labels

from ..astrometry import WCS
from ..utilities import headers
from ...logger import LOGGER

class Source(WCS):
    def __init__(self,img,seg,segid,whttype='pixels'):
        self.segid=segid
        self.whttype=whttype
        self.spectralregions=[]


        y,x=np.where((img>0) & (seg==self.segid))
        self.npixels=len(x)
        if self.npixels>0:
            # compute and subtract the local sky background:
            self.background=self.compute_local_background(img,x,y)
            img-=self.background
            
            
            # process the pixels based on the whttype specified
            if self.whttype=='pixels':
                w=np.abs(img[self.y,self.x])
            
            else:
                self.whttype='pixels':
                w=np.abs(img[self.y,self.x])

            # normalize the weights
            self.norm=np.sum(w)
            w /= self.norm
                
            # with weights process if for different segmentation regions
            # segid < 0 ---> compound
            # segid > 0 ---> simple
            #
            # Eventually have "complex" morphologies in the spectral regions
                
            if self.segid < 0:
                for xx,yy,ww in zip(x,y,w):
                    self.spectralregions.append(SpectralRegion(xx,yy,ww))
                
            elif self.segid >0:
                self.spectralregions.append(SpectralRegion(x,y,w))


            # compute the centeroids
            xyc=(np.average(x,weights=w),np.average(y,weights=w))
            self.xyc=self.image_coordinates(*xyc)
            self.adc=self.xy2ad(*xyc)

            # compute some fluxes
            self.flux=np.sum(img[y,x])
            self.fnu=self.flux*10.**(-0.4*(self.zeropoint+48.6))
            self.mag=-2.5*np.log10(self.flux)+self.zeropoint
            

    def __bool__(self):
        return bool(self.spectralregions)

    def __iter__(self):
        for region in self.spectralregions:
            yield from region

    def __getitem__(self,k):
        return self.spectralregions[k]

    def items(self):
        yield from enumerate(self.spectralregions)

    
    #def __setitem__(self,k,v):
    #    self.spectralregions[k]=v

            
            


    @property
    def name(self):
        return str(self.segid)

    #@property
    #def xd(self):
    #    xd=[]
    #    for region in self.spectralregions:
    #        xd.extend(region.x)
    #            
        


    @property
    def nregions(self):
        return len(self.spectralregions)


        

    def update_HDU(self,hdu,group=0,**kwargs):
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
        hdu.header['NPIXELS']=(self.npixels,'total number of extracted pixels')
        hdu.header['NREGIONS']=(self.nregions,'number of spectral regions')
        hdu.header['AREA']=(self.npixels*self.pixelarea,'source area (arcsec2)')
        headers.add_stanza(hdu.header,'Source Properties',before='SEGID')
        
        for k,v in kwargs.items():
            hdu.header[k]=v

        return hdu


    
    def image_coordinates(self,x,y,dtype=None):
        xx=x-self.ltv[0]
        yy=y-self.ltv[1]
        if dtype is not None:
            xx=dtype(xx)
            yy=dtype(yy)
        return xx,yy

    #def instrumental_flux(self,img):
    #    xd,yd=self.image_coordinates(self.x,self.y,dtype=np.uint32)
    #    return np.sum(img[yd,xd])


    def load_seds(self,sedlib):
        for regid,region in self.items():
            region.sed=sedlib[(segid,regid)]
            
    def write_seds(self,path='.'):
        for regid,region in self.spectralregions:
            filename=os.path.join(path,f'{self.segid}_{regid}.sed')
            region.write_sed(filename)
            
    
    def compute_local_background(self,img,x,y):
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
            return ave
        else:
            return 0.0
            
