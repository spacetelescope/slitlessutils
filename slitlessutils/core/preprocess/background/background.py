from astropy import convolution
from astropy.io import fits
from astropy.modeling import fitting,models
from astropy.stats import sigma_clip,sigma_clipped_stats
from contextlib import nullcontext
import numpy as np
import os
from scipy.signal import savgol_filter
from skimage import morphology,measure

import matplotlib.pyplot as plt

from ....logger import LOGGER
from ....config import Config
from ...utilities import headers,indices
from ...wfss import WFSS

def background_processing(mastersky=False):
    """
    A decorator so that all background methods prep the data in the same way

    Parameters
    ----------
    mastersky : bool, optional
        flag to read master sky image

    """
    
    

    
    def background(func):
        def wrapper(self,filename,newfile=None,inplace=False,**kwargs):
            


            
            # a flag on how to read the image
            if inplace:
                mode='update'
            else:
                mode='readonly'

            # load it as an observed image
            data=WFSS.observed(filename)
            if mastersky:
                # get the name of the background file:
                backfile=data.background_filename('master')

                
            with fits.open(filename,mode=mode) as hdul:
                # will need the primary header below
                phdr=hdul[0].header
                for hdu in hdul:

                    # only work on SCI images
                    if hdu.header.get('EXTNAME')=='SCI':
                        vers=hdu.header.get('EXTVER')

                        # read all the data
                        sci=hdul[('SCI',vers)].data
                        hdr=hdul[('SCI',vers)].header
                        unc=hdul[('ERR',vers)].data
                        dqa=hdul[('DQ',vers)].data

                        # get an estimate of the sky from statistics
                        gpx=(dqa==0)      # these are good pixels
                        ave,med,sig=sigma_clipped_stats(sci[gpx],
                                                        sigma=self.skysigma)

                        # if doing master sky, then need to prep the sky 
                        if mastersky:

                            # get the model and check normalization
                            mod=fits.getdata(backfile,('SKY',vers))

                            a,m,s=sigma_clipped_stats(mod)
                            if np.abs(a-1)>1e-3:
                                msg='The master sky image ({backfile}) is '+\
                                    'unnormalized. Results will be fine,\n'+\
                                    'but the units will not make sense.'
                                LOGGER.warning(msg)


                            # excise if need-be
                            if phdr['INSTRUME']=='WFC3':
                                x0=-int(hdr.get('LTV1',0))
                                y0=-int(hdr.get('LTV2',0))
                                nx=hdr['NAXIS1']
                                ny=hdr['NAXIS2']
                                mod=mod[y0:ny+y0,x0:nx+x0]

                            ret,src=func(self,sci,hdr,unc,med,gpx,mod,**kwargs)
                        else:
                            ret,src=func(self,sci,hdr,unc,med,gpx,**kwargs)
                        
                        # have the model now, so subtract it
                        hdul[('SCI',vers)].data=sci-ret
                        hdul[('SCI',vers)].header=hdr


                # do something for each type of inplace writing, but always
                # return a filename
                if inplace:
                    outfile=filename
                else:
                    if newfile is None:
                        newfile=f'{data.dataset}_sub.fits'
                    hdul.writeto(newfile,overwrite=True)
                    outfile=newfile
                                    
                
            return outfile
        return wrapper
    return background


class Background:
    """
    Class to perform background subtraction of Wide-Field slitless
    spectrosocpy.

    Parameters
    ---------
    dispaxis : str, optional
       The approximate dispersion axis, used for three purposes:

       1) set the approximate orientation of a convolution kernel that
          helps identify spectral traces
       2) set the approximate orientation of a dilation kernel that
          grows the spectral traces
       3) explicitly used by the `poly1d` method to iterate along the
          dispersion axis and fit profiles in the cross dispersion.

       The default value is 'x'.  WARNING: THIS MAY GO AWAY AS THE
       `disperser` OBJECTS KNOW ABOUT THIS AXIS.

    skysigma : float or int, optional
       The sigma clipping threshold to get an initial (assumed spatially
       flat) guess for the sky level.  Default is 3.

    srcsigma : float or int, optional
       The sigma-clipping threshold to identify spectral traces.  Default is 3

    outliersigma : float or int, optional
       The sigma-clipping threshold for fitting models in the
       `astropy.models.FittingWithOutlierRemoval()`.  Default is 3.

    minpix : int, optional
       Minimum number of pixels for a positive deviation to be considered
       significant.  Default is 5.

    maxiter : int, optional
       Maximum number of iterations for source identification    


    """

    
    MINUNC=1e-10    # minimum viable uncertainty
    
    def __init__(self,dispaxis='x',skysigma=3.,srcsigma=3.,
                 outliersigma=3.,minpix=5,maxiter=10):


        
        kernsig=0.8
        kernlen=50
        
        # make a convolution kernel
        L=4*kernsig
        x=np.linspace(-L,L,31,endpoint=True)
        kernel=np.tile(np.exp(-0.5*(x/kernsig)**2),(kernlen,1))
        kernel/=np.sum(kernel)
        

        
        self._dispaxis=dispaxis
        if self._dispaxis=='x':
            self.dispaxis=1            
            self.kernel=kernel.T
            self.footprint=morphology.rectangle(100,5)
        elif self._dispaxis=='y':
            self.dispaxis=0
            self.kernel=kernel
            self.footprint=None
        self.crossaxis=1-self.dispaxis
            
        self.skysigma=skysigma
        self.srcsigma=srcsigma
        self.outliersigma=outliersigma
        self.minpix=minpix
        self.maxiter=maxiter

    @property
    def convolve(self):
        return self.kernel is not None

    @property
    def remove_small(self):
        return self.minpix>0
    
    @property
    def dilate(self):
        return self.footprint is not None

        
    def skypixels(self,sci,unc,mod):
        """
        Simple method to find sky pixels

        Parameters
        ----------
        sci : `np.ndarray` 
           The science image to find sky pixels

        unc : `np.ndarray`
           The uncertainty image for sigma-thresholding

        mod : float or `nd.ndarray`
           The model sky image.  Must be broadcastable to the shape
           of `sci` and `unc`.

        Returns
        -------
        sky : `np.ndarray` (bool dtype)
           The pixels that belong to the sky
        """
        

        # convolve?
        if self.convolve:
            con=convolution.convolve_fft(sci,self.kernel)
        else:
            con=sci
    
        # sigma thresholding for sources (ie. pix brighter than some limit)
        src=(con-mod)>(self.srcsigma*unc)

        # remove small objects
        if self.remove_small:
            src=morphology.remove_small_objects(src,min_size=self.minpix)
                    
        # dilate
        if self.dilate:
            src=morphology.dilation(src,footprint=self.footprint)

        # return the sky pixels (ie those that are *NOT* sources)
        return np.logical_not(src)

    @background_processing(mastersky=True)
    def master(self,sci,hdr,unc,mod,gpx,img):
        """
        Function to fit a *single* master sky image to a WFSS image.

        The master-sky image is a two-dimensional model of the dispersed
        sky background.  Ideally, this image is scaled to match the background
        pixels in a dispersed image, however requires flagging all of the 
        pixels that contain any additional signal (such as, but not limited to, 
        astrophysical sources from any spectral order, cosmic rays, 
        persistence or any type of bad pixels).  The masking is implemented
        via an iterative approach:
        
        0) initalize the background model as the 3-sigma clipped median 
           (:math:`m`) of the pixels not flagged in the data-quality array 
            and flag pixels that are 
        
        .. math::  
           \\frac{{\cal I}-m}{{\cal U}} \leq n_{sig}
    
        where :math:`{\cal I}` and :math:`{\cal U}` are the science and 
        uncertainty images, respectively.    

        1) compute optimal scaling coefficient as

        .. math::
           \\alpha = \\frac{F_{OT}}{F_{TT}}

        where 

        .. math::
           :nowrap:
           \begin{eqnarray}
             F_{OO} & = & \sum_{x,y\in s} \left(\\frac{{\cal I}}{{\cal U}}\\right)^2\\
             F_{OT} & = & \sum_{x,y\in s} \\frac{{\cal I}}{{\cal U}}\\frac{{\cal S}}{{\cal U}}\\
             F_{TT} & = & \sum_{x,y\in s} \left(\\frac{{\cal S}}{{\cal U}}\\right)^2\\
           \end{eqnarray}
    
        and :math:`s` is the collection of pixels that are neither flagged in 
        the data-quality arrays nor above the background level.

        2) redefine the non-background pixels by testing:
        .. math:: 
           \\frac{{\cal I}-\\alpha{\cal S}}{{\cal U}} \leq n_{sig}

        3) go to step 1) until either the maximum number of iterations is
        reached or the fractional difference is less than the set tolerance 
        (:math:`\epsilon`):

        .. math::
          |\\alpha_{i-1} - \\alpha_i| \leq \epsilon |\\alpha_{i}|

        Parameters
        ----------
        data : str or `wfss.WFSS` 
           Either a string that is the full path to a WFSS file or a 
           `wfss.WFSS`.  


        Returns
        -------
        outfile : str
           The name of the master-sky subtracted file.
    
        
        """

        
        # do inverse-variance weighting:
        sci2=sci/np.maximum(unc,self.MINUNC)
        img2=img/np.maximum(unc,self.MINUNC)

        # find the sky pixels
        sky=self.skypixels(sci,unc,mod)
        
        # start iteration
        npix=np.count_nonzero(sky)
        proceed=True
        it=1
        while proceed:

            # find the good pixels
            g=np.logical_and(sky,gpx)
            s2=sci2[g]
            i2=img2[g]
            
            # compute some auxillary variables for linear fitting
            Foo=np.sum(s2*s2)
            Fot=np.sum(s2*i2)
            Ftt=np.sum(i2*i2)
            
            # the scale facto
            alpha=Fot/Ftt

            # the chi2 and best model
            chi2=Foo-alpha*Fot
            out=alpha*img

            # update the sky pixel mask            
            sky=self.skypixels(sci,unc,out)

            # update iteration variables
            it+=1
            npix2=np.count_nonzero(sky)
            proceed=(npix != npix2) and (it < self.maxiter)
            npix=npix2
            if not proceed:
                break

        if alpha < 0:
            msg=('The expected sky is negative. This is usually because:',
                 '  1) the image is already background subtracted, or',
                 '  2) there is some catastrophic anamoly or large source')
            LOGGER.warning('\n'.join(msg))
                

            
        # update the header
        self.update_header(hdr)
        hdr.set('METHOD',value='master',comment='method')
        hdr.set('ALPHA',value=alpha,
                comment='Normalization of sky model')
            
            
        if it> self.maxiter:
            LOGGER.warning('maximum number iterations exceeded.')

        return out,np.logical_not(sky)

         
    

    @background_processing(mastersky=False)
    def poly1d(self,sci,hdr,unc,mod,gpx,degree=2,filtwindow=51,filtorder=1):

        """
        Method to fit polynomials in the cross-dispersion axis and
        smooth with Savitzky-Golay in the dispersion axis.

        The algorithm iterates over the dispersion axis and fits polynomials
        in the cross dispersion axis.  The fitting is done with
        `astropy.models.FittingWithOutlierRejection`, which creates a
        temporary background model.  Since this model will have significant
        pixel-to-pixel noise, a perpendicular axis is iterated over with
        a Savtizky-Golay filter.

        Parameters
        ----------
        data : str or `wfss.WFSS` 
           Either a string that is the full path to a WFSS file or a 
           `wfss.WFSS`.  

        degree : int, optional
           Order of the cross-dispersion polynomial.  Default is 2

        filtwindow : int, optional
           The window size of the Savitzky-Golay filter.  Default is 51

        filtorder : int, optional
           The polynomial order of the Savitzky-Golay filter, which is
           applied in the disperion axis.  Default is 1
        
        Returns
        -------
        outfile : str
           The name of the master-sky subtracted file.
    

        Notes
        -----
        Below we assume that :math:`\lambda` and :math:`\eta` are the
        dispersion and cross-dispersion axes, respectively.  
        """
        
        sky=self.skypixels(sci,unc,mod)
        fitter=fitting.LinearLSQFitter()
        ofitter=fitting.FittingWithOutlierRemoval(fitter,sigma_clip,
                                                  sigma=self.outliersigma)

        # the output model
        out=np.zeros_like(sci,dtype=float)

     
        # set some slice objects
        if self.dispaxis==1:
            loopdisp=lambda lam: (slice(None,None,None),lam)
            loopcross=lambda eta: (eta,slice(None,None,None))
        elif self.dispaxis==0:
            loopcross=lambda lam: (slice(None,None,None),lam)
            loopdisp=lambda eta: (eta,slice(None,None,None))
            
        

        # set up the iteration
        npix=np.count_nonzero(sky)
        proceed=True
        it=1
        while proceed:

            # get estimate by fitting in cross dispersion axis
            eta=np.linspace(-1,1,sci.shape[self.crossaxis],dtype=float)
            model=models.Polynomial1D(degree)
            for lam in range(sci.shape[self.dispaxis]):
                # get pixels in the cross dispersion axis
                xy=loopdisp(lam)                

                # compute weights and do a fit
                w=np.logical_and(sky[xy],gpx[xy])/unc[xy]**2
                pars,mask=ofitter(model,eta,sci[xy],weights=w)

                # update the model
                out[xy]=pars(eta)


            # smooth estimates in the dispersion axis
            for eta in range(sci.shape[self.crossaxis]):
                xy=loopcross(eta)
                out[xy]=savgol_filter(out[xy],filtwindow,filtorder,
                                      deriv=0,mode='nearest')
                                      
                                         
                
            # update sky pixels
            sky=self.skypixels(sci,unc,out)

            # update iteration variables
            it+=1
            npix2=np.count_nonzero(sky)
            proceed=(npix > npix2) and (it < self.maxiter)
            npix=npix2


        # update the header
        self.update_header(hdr)
        hdr.set('METHOD',value='poly1d',comment='method')
        hdr.set('CROSSORD',value=degree,
                comment='Polynomial order in cross dispersion')
        hdr.set('FILTWIND',value=filtwindow,
                comment='Savitzky-Golay filter window size')
        hdr.set('FILTORD',value=filtorder,
                comment='Savitzky-Golay filter order')
        hdr.set('ITER',value=it,comment='number of iterations')
        hdr.set('NPIX',value=npix,comment='number of sky pixels')

        

            
        if it> self.maxiter:
            LOGGER.warning('maximum number iterations exceeded.')
            
        return out,np.logical_not(sky)

    def update_header(self,hdr):
        """
        Method to update an `astropy.fits.Header()` object

        Parameters
        ---------
        hdr : `astropy.fits.Header`
            The fits header to update
        """

        
        hdr.set('SKYNSIG',value=self.skysigma,
                comment='Nsigma for constant sky model')

        hdr.set('SRCNSIG',value=self.srcsigma,
                comment='Nsigma for thresholding sources')
        hdr.set('OUTNSIG',value=self.outliersigma,
                comment='Nsigma for astropy.model outlier rejection')
        hdr.set('MINPIX',value=self.minpix,
                comment='minimum number of pixels to be source')
        hdr.set('MAXITER',value=self.maxiter,
                comment='Maximum number of iterations for source flagging')
        headers.add_stanza(hdr,'Background Subtraction',before='SKYNSIG')
        
if __name__=='__main__':
    #plt.ion()
    #fig = plt.figure()
    #ax = fig.add_subplot(111)
    #ax.imshow(out)
    #fig.canvas.draw()
    #fig.canvas.flush_events()
    #plt.pause(0.1)
    pass
