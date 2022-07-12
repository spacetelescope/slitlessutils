import numpy as np
import os
from skimage.segmentation import expand_labels
import matplotlib.pyplot as plt
from astropy.stats import sigma_clipped_stats


from .spectralregion import SpectralRegion
from ..astrometry import WCS
from ..utilities import headers,indices
from ...logger import LOGGER

class Source(WCS):
    def __init__(self,img,seg,segid,local_back=True,backsize=5,nsig=(5,5),
                 whttype='pixels',zeropoint=25.):
        WCS.__init__(self,img.header)

        # record some things
        self.segid=segid
        self.whttype=whttype
        self.backsize=backsize
        self.local_back=local_back

        # make an empty list for regions
        self.spectralregions=[]

        # test the number of sigma for the sky
        if isinstance(nsig,(float,int)):
            nsig=(nsig,nsig)
        self.nsig=nsig



        # find the pixels for this segmentation
        y,x=np.where((img>0) & (seg==self.segid))
        self.npixels=len(x)


        if self.npixels>0:
            # compute and subtract the local sky background:
            if self.local_back and self.backsize>0:
                # get a bounding box
                x0=max(np.amin(x)-self.backsize,0)
                y0=max(np.amin(y)-self.backsize,0)
                x1=min(np.amax(x)+self.backsize,img.shape[1]-1)
                y1=min(np.amax(y)+self.backsize,img.shape[0]-1)

                # cut out regions
                subimg=img.image[y0:y1,x0:x1]
                subseg=seg.image[y0:y1,x0:x1]

                # grow the segmapsize the segmap
                subseg2=expand_labels(subseg,distance=self.backsize)
                g=np.where((subseg2==segid) & (subseg==0))

                # compute average in background region
                ave,med,sig=sigma_clipped_stats(subimg[g],
                                                sigma_lower=self.nsig[0],
                                                sigma_upper=self.nsig[1])
                #ave=np.average(subimg[g])
                self.background=ave
            else:
                self.background=0.0

            img-=self.background


            # process the pixels based on the whttype specified
            if self.whttype=='pixels':
                w=np.abs(img.image[y,x])

            else:
                self.whttype='pixels'
                w=np.abs(img.image[y,x])

            # normalize the weights
            self.norm=np.sum(w)
            w /= self.norm


            # put coordinates back on original footprint
            #x,y=self.image_coordinates(x,y,dtype=np.int)


            #
            # with weights process if for different segmentation regions
            # segid < 0 ---> compound
            # segid > 0 ---> simple
            #
            # Eventually have "complex" morphologies in the spectral regions
            # that can describe an arbitrary arrangement of spectral regions
            #
            if self.segid < 0:
                # a compound source, so each pixel in a separate SpectralRegion
                for xx,yy,ww in zip(x,y,w):
                    self.spectralregions.append(SpectralRegion([xx],[yy],[ww]))
            elif self.segid >0:
                self.spectralregions.append(SpectralRegion(x,y,w))


            # compute the centeroids
            xyc=(np.average(x,weights=w),np.average(y,weights=w))
            self.xyc=self.image_coordinates(*xyc)
            self.adc=self.xy2ad(*xyc)


            # compute the area of the object
            self.area=self.npixels*self.pixelarea



            # compute some fluxes
            self.flux=np.sum(img.image[y,x])
            self.fnu=self.flux*10.**(-0.4*(zeropoint+48.6))
            self.mag=-2.5*np.log10(self.flux)+zeropoint


    def __bool__(self):
        return bool(self.spectralregions)

    def __getitem__(self,k):
        return self.spectralregions[k]

    def __iter__(self):
        for region in self.spectralregions:
            yield from region

    def __next__(self):
        if not hasattr(self,'_iter'):
            self._iter=iter(self)
        return next(self._iter)

    def __str__(self):
        return super().__str__()


    #def image_pixels(self,dtype=None):
    #    for region in self.spectralregions:
    #        x,y=source.image_coordiantes(region.x,region.y,dtype=dtype)
    #        yield from zip(x,y)

    #def source_pixels(self):
    def pixels(self):
        for region in self.spectralregions:
            yield from region.pixels()

    def items(self):
        yield from enumerate(self.spectralregions)

    @property
    def name(self):
        return str(self.segid)

    @property
    def nregions(self):
        return len(self.spectralregions)


    def update_header(self,hdr,group=0):

        hdr['SEGID']=(self.segid,'Segmentation ID')
        hdr['GRPID']=(group,'Group ID')
        hdr['RA']=(self.adc[0],'Right Ascension (deg)')
        hdr['DEC']=(self.adc[1],'Declination (deg)')
        hdr['X']=(self.xyc[0],'X barycenter')
        hdr['Y']=(self.xyc[1],'Y barycenter')
        hdr['FLUX']=(self.flux,'instrumental flux in direct image')
        hdr['MAG']=(self.mag,'AB magnitude in direct image')
        hdr['FNU']=(self.fnu,'flux in erg/s/cm2/Hz in direct image')
        hdr['NPIXELS']=(self.npixels,'total number of extracted pixels')
        hdr['WHTTYPE']=(self.whttype,'Type of source profile weights')
        hdr['NREGIONS']=(self.nregions,'number of spectral regions')
        hdr['AREA']=(self.npixels*self.pixelarea,'source area (arcsec2)')
        headers.add_stanza(hdr,'Source Properties',before='SEGID')


        hdr['BCKSUB']=(self.local_back,'Was local background in direct image subtracted')
        hdr['BCKSIZE']=(self.backsize,'Size of background annulus (in pix)')
        hdr['BCKVAL']=(self.background,'Background level in direct image')
        hdr['BCKLOSIG']=(self.nsig[0],'N sigma for low level')
        hdr['BCKHISIG']=(self.nsig[1],'N sigma for high level')
        headers.add_stanza(hdr,'Direct Image Background',before='BCKSUB')



        #for k,v in kwargs.items():
        #    hdr[k]=v



    def image_coordinates(self,x,y,dtype=None):
        xx=x-self.ltv[0]
        yy=y-self.ltv[1]
        if dtype is not None:
            xx=xx.astype(dtype)
            yy=yy.astype(dtype)

            #xx=dtype(xx)
            #yy=dtype(yy)
        return xx,yy

    #def instrumental_flux(self,img):
    #    ''' this doesnt make a lot of sense for ground-based images '''
    #    xd,yd=self.image_coordinates(self.x,self.y,dtype=np.uint32)
    #    return np.sum(img[yd,xd])


    def load_seds(self,sedlib,throughput=None):
        for regid,region in self.items():
            sed=sedlib[(self.segid,regid)]
            if throughput:
                sed.normalize(throughput,self.fnu*np.sum(region.w))

            region.sed=sed

    def write_seds(self,path='.',**kwargs):
        for regid,region in enumerate(self.spectralregions):
            filename=os.path.join(path,f'{self.segid}_{regid}.sed')
            region.sed.write_file(filename,**kwargs)




    def plot(self):

        x=[]
        y=[]
        w=[]
        r=[]
        for regid,region in enumerate(self.spectralregions):
            x.extend(region.x)
            y.extend(region.y)
            w.extend(region.w)
            r.extend([regid]*len(region))

        x=np.array(x)
        y=np.array(y)
        w=np.array(w)
        r=np.array(r)


        x0,x1=np.amin(x),np.amax(x)
        y0,y1=np.amin(y),np.amax(y)
        shape=(y1-y0+1,x1-x0+1)

        wht=np.full(shape,np.nan)
        reg=np.full(shape,np.nan)
        wht[y,x]=w
        reg[y,x]=r


        #
        #
        # x0,x1=np.inf,0
        # y0,y1=np.inf,0
        #
        # for region in self.spectralregions:
        #     x0=min(x0,np.amin(region.x))
        #     x1=max(x1,np.amax(region.x))
        #     y0=min(y0,np.amin(region.y))
        #     y1=max(y1,np.amax(region.y))
        #
        # shape=(y1-y0+1,x1-x0+1)
        # wht=np.full(shape,np.nan,dtype=float)
        # reg=np.full(shape,np.nan,dtype=float)
        #
        # for regid,region in enumerate(self.spectralregions):
        #     wht[region.y,region.x]=region.w
        #     reg[region.y,region.x]=regid


        fig,axes=plt.subplots(1,2,sharex=True,sharey=True)
        axes[0].imshow(wht,origin='lower',interpolation='nearest')
        axes[0].set_title('weights')
        axes[1].imshow(reg,origin='lower',interpolation='nearest')
        axes[1].set_title('region')
        plt.tight_layout()
        fig.canvas.draw()
        plt.show()
