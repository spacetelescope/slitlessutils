import slitlessutils as su
import os
import pandas as pd
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.modeling import models
from skimage.morphology import label


from .parameters import *




def mkscene(**kwargs):
  
    # get the SED path
    if 'PYSYN_CDBS' in os.environ:
        sedpath=os.path.join(os.environ['PYSYN_CDBS'],'grid',
                             'pickles','dat_uvk')
        if not os.path.isdir(sedpath):
            su.LOGGER.error("Required: Pickles SED library")
            return
    else:
        su.LOGGER.error('Required: the PYSYN_CDBS set')
        return

    dataset=['_'.join((dset,ROOT)) for dset in ('a','b','c')]
    ndataset=len(dataset)
    orientat=np.round(np.arange(0,180,180./ndataset),0)
    
    # establish the pointing of the simulated grism images
    df=pd.DataFrame(data={'dataset':dataset,
                          'ra': [RA]*ndataset,
                          'dec':[DEC]*ndataset,
                          'orientat':orientat,
                          'telescope':['HST']*ndataset,
                          #'instrument':['WFC3IR','WFC3IR','WFC3IR'],
                          #'grating':['G102','G102','G102'],
                          'instrument':['ACSWFC']*ndataset,
                          'grating':['G800L']*ndataset,
                          'pupil':['']*ndataset})
    df.set_index('dataset',inplace=True,drop=True)
    df.to_csv(f'{ROOT}_wcs.csv')

    
    # make direct images.  start by filling a direct image
    img=np.zeros((NPIX,NPIX),dtype=float)
    xim,yim = np.meshgrid(np.arange(NPIX), np.arange(NPIX))
    
    # make the direct image morphology the same for all sources:
    func = models.Gaussian2D(x_stddev=1.5,y_stddev=1.5)

    # properties of the sources
    x=[500,520,420]   # x-centers
    y=[500,500,550]   # y-centers
    a=[2.,1.,0.5]     # amplitudes (ie. peak fluxes)
    f=['pickles_uk_21','pickles_uk_31','pickles_uk_41']


    x=[500]   # x-centers
    y=[500]   # y-centers
    a=[2.,]     # amplitudes (ie. peak fluxes)
    f=['pickles_uk_21']


    
    # create an empty SED HDUList
    hdul=fits.HDUList()
    
    # put the sources in the image and update SEDs
    # NB: start=1 here is because segid=0 is reserved for sky
    for i,(x_mean,y_mean,amp,sfile) in enumerate(zip(x,y,a,f),start=1):
        func.x_mean=x_mean
        func.y_mean=y_mean
        func.amplitude=amp
        img+= func(xim,yim)

        # create empty header for SED and add it to the SEDfile
        hdr=fits.Header()
        hdr['EXTNAME']=(str(i),'segmentation ID')
        hdr['EXTVER']=(0,'spectral region ID')
        hdr['FILENAME']=os.path.join(sedpath,sfile)+'.fits'
        hdul.append(fits.BinTableHDU(header=hdr))

        
    # make segmentation map
    threshold = 0.2                # threshold to apply to the image
    seg = label(img > threshold)   # label the sources with unique seg IDs
            
    # make a header:
    w=WCS(naxis=2)
    w.wcs.crpix=[NPIX/2.,NPIX/2.]
    w.wcs.crval=[RA,DEC]
    w.wcs.ctype=['RA---TAN','DEC--TAN']
    w.wcs.cd=[[-PIXSCL/3600.,0.],[0.,PIXSCL/3600.]]
    hdr=w.to_header()


    # add some info on the filter
    #hdr['FILTFILE']='hst_acs_f775w.fits'
    hdr['TELESCOP']='hst'
    hdr['INSTRUME']='acs'
    hdr['FILTER']='f775w'

    
    # save to disk
    fits.writeto(f'{ROOT}_seg.fits',seg,header=hdr,overwrite=True)
    fits.writeto(f'{ROOT}_img.fits',img,header=hdr,overwrite=True)
    hdul.writeto(f'{ROOT}_seds.fits',overwrite=True)

