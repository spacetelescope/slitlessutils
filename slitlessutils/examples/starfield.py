import slitlessutils as su

from astropy.io import fits
from astropy.wcs import WCS
from astropy.modeling import models
from skimage.morphology import label
import matplotlib.pyplot as plt
import numpy as np

import os

from .parameters import *


ROOT='starfield'            # base name for files in this example
DATASETS=('a','b','c')      # dataset names (will do 3 datasets)
TELESCOPE='HST'             # name of the telescope
INSTRUMENT='ACS'            # name of the instrument
DETECTOR='WFC'              # name of the detector in the instrument
DISPERSER='G800L'           # name of the spectral element
BLOCKING=''                 # name of the blocking filter
WRANGE=(5500.,9500.)        # wavelength range to inspect
FILTER='F775W'              # name of the filter for the direct image
SUFFIX = 'flc'              # type of suffix
NCPU=1                      # number of CPUs to use


DETECTOR='SBC'
DISPERSER='P130L'
WRANGE=(1150.,1900.)
FILTER='F125LP'
SUFFIX='flt'






SPECCAT='pickles'

# properties of the sources.  The entries are: SEGID:(x,y,ABmag,sedfile)
SOURCES={1:(500,500,19.0,'pickles_uk_15.fits'),
         2:(485,480,19.5,'pickles_uk_20.fits'),
         3:(420,550,20.0,'pickles_uk_25.fits')}
PSFSIG=1.5          # stdev of the PSF in pixels
APERRAD=0.3         # radius in arcsec     
SEDPATH='starfield_input'


"""
A set of routines to demonstrate the creation and analysis of WFSS for a 
simple scene of point sources, each having a unique spectrum and normalization.


"""




def make_scene():
    """
    Function to make a scene of point sources.

    Notes
    -----
    1) This is to make direct image and segmentation map, as such, this is 
       not likely needed to demonstrate slitlessutils, since these would 
       likely be synthesized by other means (such as astrodrizzle, source 
       extractor, and/or photutils).  Therefore, this is just to simulate 
       those products.

    2) Assumes a Gaussian for the point-source function (PSF).  This is not 
       a critical assumption, but rather one done for simplicity.

    3) Writes all products to disk.

    """
            
    # create something to store the SEDs that will be used for the simulation
    hdul=fits.HDUList()

    # make a direct and segmentation image
    img=np.zeros((NPIX,NPIX),dtype=float)
    seg=np.zeros_like(img,dtype=int)
    xim,yim=np.meshgrid(np.arange(NPIX),np.arange(NPIX))

    # create a morphological function for the sources
    func=models.Gaussian2D(x_stddev=PSFSIG,y_stddev=PSFSIG)


    # read the filter curve
    band=su.photometry.Throughput.from_keys(TELESCOPE,INSTRUMENT,FILTER)
    
    # populate the image for each source
    for segid,(x,y,m,f) in SOURCES.items():

        # download the file from CDBS
        localfile=su.photometry.SED.get_from_CDBS(SPECCAT,f)

        
        # put the spectrum in the fits file
        hdr=fits.Header()
        hdr['EXTNAME']=(str(segid),'segmentation ID')
        hdr['EXTVER']=(0,'spectral region')
        hdr['FILENAME']=f
        hdul.append(fits.BinTableHDU(header=hdr))
        
                
        # problem here.  Astropy uses "amplitude" as the overall
        # normalization of the spatial profile, even though what we want
        # is the total flux to be specified.  Of course for some analytic
        # profiles, this amplitude can be determined, however I want
        # to leave the profile flexible, so will determine the conversion
        # between amplitude and total by just scaling via a "delta"
        # image (called dim)

        # evaluate the source profile
        func.x_mean=x
        func.y_mean=y
        func.amplitude=1.  # set the amplitude to 1, renormalize later
        dim=func(xim,yim)
        
        # get the expected total flux given the requested AB mag
        ftot=10.**(-0.4*(m-band.zeropoint))

        # scale the delta image and add into the image
        img+=(dim*(ftot/np.sum(dim)))
        

        # compute the aperture flux
        rim=np.hypot(xim-x,yim-y)
        g=np.where(rim < APERRAD/PIXSCL)
        seg[g]=segid
        
        
        
    # make header for the direct and segmentation image
    w=WCS(naxis=2)
    w.wcs.crpix=[NPIX/2.,NPIX/2.]
    w.wcs.crval=[RA,DEC]
    w.wcs.ctype=['RA---TAN','DEC--TAN']
    w.wcs.cd=[[-PIXSCL/3600.,0.],[0.,PIXSCL/3600.]]
    h=w.to_header()

    # add some info to the file
    h['TELESCOP']=TELESCOPE
    h['INSTRUME']=INSTRUMENT
    h['FILTER']=FILTER
    h['ZERO']=band.zeropoint
    
    # write the images to disk
    fits.writeto(f'{ROOT}_seg.fits',seg,header=h,overwrite=True)
    fits.writeto(f'{ROOT}_sci.fits',img,header=h,overwrite=True)
    hdul.writeto(f'{ROOT}_seds.fits',overwrite=True)
    


def simulate_grisms():
    """
    Method to simulate the grisms

    Notes
    -----
    1) This is the primary demonstration of how to simulate grism image(s).   
    
    2) The notional observational setup will be specified in a csv file. 
       The one needed here is made here, but shows how one can make their 
       own for some other purpose.
    
    3) All products are written to disk.
    
    """


    # write a "WCS" file to disk that contains properites of
    # the images to emulate an observers setup
    with open(f'{ROOT}_wcs.csv','w') as fp:
        print('dataset,ra,dec,orientat,telescope,instrument,disperser,blocking',file=fp)

        n=len(DATASETS)

        # put each WFSS exposure in the canon
        for i,dset in enumerate(DATASETS):
            orientat=i*(180./n)
            line=','.join((dset+'_'+ROOT,str(RA),str(DEC),str(orientat),
                           TELESCOPE,INSTRUMENT+DETECTOR,DISPERSER,BLOCKING))
            print(line,file=fp)

    
    # load the grism images
    data=su.wfss.WFSSCollection.from_wcsfile(f'{ROOT}_wcs.csv')
    
    # load the sources
    sources=su.sources.SourceCollection(f'{ROOT}_seg.fits',f'{ROOT}_sci.fits',
                                        sedfile=f'{ROOT}_seds.fits')
    
    # save the normalized SEDs
    sources.write_seds(path=SEDPATH)

    
    # project the sources onto the grism images
    tab=su.modules.Tabulate(ncpu=NCPU)
    pdtfiles=tab(data,sources)

    # use the projection tables and the SEDs to simulate grism images
    sim=su.modules.Simulate(ncpu=NCPU)
    imgfiles=sim(data,sources)

    return imgfiles


    
def extract_single():
    """
    Method to demonstrate single-exposure spectral extraction

    Notes
    -----
    1) Normally, one would need to "tabulate" the images before extracting,
       but this should've been done after the simulations.  If this needs
       to be done, then see `simulate_grisms()` for an example.
    
    2) Will write a file to disk with name: {ROOT}_x1d.fits

    """

    # load the grism images
    data=su.wfss.data.WFSSCollection.from_glob(f'*{ROOT}_{SUFFIX}.fits.gz')

    # load the sources into SU
    sources=su.sources.SourceCollection(f'{ROOT}_seg.fits',
                                        f'{ROOT}_sci.fits')

    
    # run the single-file extraction
    ext = su.modules.Single('+1',root=ROOT,ncpu=NCPU)
    res = ext(data,sources)
    return res


def extract_multi():
    """
    Method to demonstrate multi-exposure spectral extraction

    Notes
    -----
    1) Normally, one would need to "tabulate" the images before extracting,
       but this should've been done after the simulations.  If this needs
       to be done, then see `simulate_grisms()` for an example.
    
    2) Will write a file to disk with name: {ROOT}_multi_x1d.fits

    """

    
    # load the grism images
    data=su.wfss.data.WFSSCollection.from_glob(f'*{ROOT}_{SUFFIX}.fits.gz')

    # load the sources into SU
    sources=su.sources.SourceCollection(f'{ROOT}_seg.fits',
                                        f'{ROOT}_sci.fits')
    

    # run the multi-orient extraction
    ext = su.modules.Multi('+1',(-3.,1.,0.1),algorithm='grid')
    res = ext(data,sources,root=ROOT+'_multi')

    return res


def extract_group():
    """
    Method to demonstrate grouping methods used in the mult-exposure 
    extraction

    Notes
    -----
    1) Normally, one would need to "tabulate" the images before extracting,
       but this should've been done after the simulations.  If this needs
       to be done, then see `simulate_grisms()` for an example.
    
    2) Will write a file to disk with name: {ROOT}_group_x1d.fits
    
    3) This method is largely the same as `extract_multi()` with some 
       additional bits for the grouping
    
    """

    # load the grism images
    data=su.wfss.data.WFSSCollection.from_glob(f'*{ROOT}_{SUFFIX}.fits.gz')

    # load the sources into SU
    sources=su.sources.SourceCollection(f'{ROOT}_seg.fits',
                                        f'{ROOT}_sci.fits')
    

    # run the grouping algorithms
    grp=su.modules.Group(orders=('+1',),ncpu=NCPU)
    groups=grp(data,sources)

    # pass the groups to the multi-extract
    ext = su.modules.Multi('+1',(-3.,1.,0.1),algorithm='grid')
    res = ext(data,sources,root=ROOT+'_group',groups=groups)


def regions():
    """
    Method to demonstrate the region-creation tool

    Notes
    -----
    1) Normally, one would need to "tabulate" the images before extracting,
       but this should've been done after the simulations.  If this needs
       to be done, then see `simulate_grisms()` for an example.
    
    2) Will write one file for each dataset and detector, they will be 
       named: '{dataset}_{detector}.reg'
        
    """
    
        
    # load data into SU
    data=su.wfss.data.WFSSCollection.from_glob(f'*{ROOT}_{SUFFIX}.fits.gz')
    
    # load the sources into SU
    sources=su.sources.SourceCollection(f'{ROOT}_seg.fits',
                                        f'{ROOT}_sci.fits')

    reg=su.modules.Region(ncpu=1)
    res=reg(data,sources)

    
    


    
    
def compare(nsig=1.):
    """
    Method to compare the inputs to the outputs

    Parameters
    ----------
    nsig : float or int, optional
        The number of sigma in the uncertainties to plot.  Default is 1.

    """

    n=len(SOURCES)
    m=2*n
    aspect=4./3.
    fig,axes=plt.subplots(n,1,figsize=(m,m*aspect),sharex=True)


    
    regid=0
    for ax,(segid,(x,y,mag,sedfile)) in zip(axes,SOURCES.items()):
          
        filename=os.path.join(SEDPATH,f'{segid}_{regid}.sed')
        l,f=np.loadtxt(filename,usecols=(0,1),unpack=True)
        f/=1e-17

        g=np.where((WRANGE[0] <= l) & ( l <= WRANGE[1]))

        ylim=(np.amin(f[g])*0.9,np.amax(f[g])*1.1)
        ax.set_ylim(*ylim)
        ax.set_xlim(*WRANGE)
        model=ax.plot(l,f,label='input spectrum',color='black')
    
    
        

    files={'single-exposure':(f'{ROOT}_x1d.fits','green'),
           'multi-exposure':(f'{ROOT}_multi_x1d.fits','blue'),
           'group-exposure':(f'{ROOT}_group_x1d.fits','red')}

    artists=[model[0]]
    labels=['model spectrum']
    for label,(f,c) in files.items():
        if os.path.exists(f):

            first=True     # flag to append plots
            with fits.open(f,mode='readonly') as hdul:                
                for ax,(segid,(x,y,mag,sedfile)) in zip(axes,SOURCES.items()):
                    s=str(segid)
                    if s in hdul:
                        hdu=hdul[s]


                        
                        
                        lo=hdu.data['flam']-nsig*hdu.data['func']
                        hi=hdu.data['flam']+nsig*hdu.data['func']
                       
                        patch=ax.fill_between(hdu.data['lamb'],lo/10000,hi/10000,
                                              color=c,alpha=0.2)
                        meas=ax.plot(hdu.data['lamb'],hdu.data['flam']/10000,color=c)

                        if first:
                            artists.append((patch,meas[0]))
                            labels.append(label)
                            first=False
                            
    axes[0].legend(artists,labels)
    axes[n//2].set_ylabel('$f_{\lambda}$ (10$^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$)')
    axes[-1].set_xlabel('wavelength ($\mathrm{\AA}$)')


    plt.tight_layout()
    plt.savefig(f'{ROOT}_comparisons.pdf')
    plt.show()


    

def run_all():
    """
    Method to run all of the substeps in a single call

    Notes
    -----
    1) See individual functions for their descriptions
    2) Order of steps:
       make_scene(), simulate_grisms(), extract_single(), extract_multi(),
       extract_group()

    """
    
    make_scene()

    simulate_grisms()

    extract_single()
    extract_multi()
    extract_group()
    #compare('multi')
    compare()
