
from astropy.io import fits
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

import os
import pwd


from .....config import Config,SUFFIXES
from ...group import GroupCollection
from .....info import __code__
from .matrix import Matrix
from ...module import Module
from .....logger import LOGGER
from .optimizer import optimizer
from ....utilities import as_iterable,headers


class Multi(Module):
    """
    Module to implement multi-ended spectral extraction, as detailed in 
    Ryan, Casertano, & Pirzkal (2018)

    inherits from `su.core.modules.Module`

    Parameters
    ----------
    extorders : str, or iterable
       The name of the spectral orders to subtract.  This is likely a scalar
       string, but could be an iterable of strings to extract multiple
       spectral orders.

    logdamp : int, float, list, tuple 
       This gets passed to the optimizer, so its properties are described
       in those routines
    
    mskorders : many types, optional
       The orders to mask when doing the mult-ended extraction, and this 
       can take on different types.  If it is a list, set, or tuple, then 
       it is interpreted as multiple orders.  If it is the string 'all', 
       then all of the orders in the configuration are used.  If it is 
       set as `None`, then no orders are masked.      

    algorithm : str, optional
       The optimization algorithm.  Default is 'golden'


    Notes
    -----
    See the paper Ryan, Casertano, & Pirzkal (2018)
    https://ui.adsabs.harvard.edu/abs/2018PASP..130c4501R/abstract
    for more details on the algorithm.  But in brief, the flux in the 
    WFSS pixels is modeled as a weighted sum over all sources, all 
    wavelengths, and all spectral orders.  Therefore if the weighting 
    elements can be determined, then the spectra can be directly 
    inferred by inverting the linear system of equations.  In practice,
    these weights are estimated in the `tabulate()` procedure, and 
    the equations are cast as a linear operator `Matrix`, which is 
    iteratively inverted using sparse-linear algebra techniques.      

    """
        
    DESCRIPTION = "Extracting (multi)"

    def __init__(self,extorders,logdamp,mskorders=None,algorithm='golden',
                 **kwargs):
        Module.__init__(self,self.extract,**kwargs,multiprocess=False)

        self.extorders=as_iterable(extorders)
        self.mskorders=as_iterable(mskorders)
        
        self.optimizer=optimizer(algorithm,logdamp)
        self.matrix=Matrix(self.extorders,**kwargs)


    def extract(self,data,sources,groups=None,root=None):
        """
        Method to do the mult-ended spectral extraction

        Parameters
        ----------
        data : `su.core.wfss.WFSSCollection`
            The collection of WFSS images

        sources : `su.core.sources.SourceCollection`
            The collection of sources
        
        groups : `GroupCollection` or None
            The collection of groups.  If set as None, then no grouping
            will be performed.  Default is None

        root : str, optional
            The root name of the output products.  If None, then 
            'slitlessutils' is used.  Default is None

        Notes
        -----
        This is likely not to be directly called.

        """

        
        # outputing structure
        hdul1=fits.HDUList()     # simple extractions
        hdul3=fits.HDUList()     # compound extractions
        hdulL=fits.HDUList()     # L-curve data

        # make a primary header
        phdu=fits.PrimaryHDU()
        headers.add_preamble(phdu.header,
                             filetype=(' ','contents of this file'),
                             exttype=('multi','method for extraction'))
                             
        headers.add_software_log(phdu.header)              # about software
        data.update_header(phdu.header)                # about WFSS images
        sources.update_header(phdu.header)                 # about the sources
        self.optimizer.update_header(phdu.header)
        Config().update_header(phdu.header)                # global config info

        # add the primary to the file
        hdul1.append(phdu)
        hdul3.append(phdu)

        # file root name
        if not isinstance(root,str):
            root=__code__
            

        # put loops over groups here
        if groups:
            LOGGER.info('Using a group catalog')
        else:
            LOGGER.info('No grouping')
            groups=GroupCollection()



        # open a PDF to write Grouping images
        pdffile=f'{root}_{SUFFIXES["L-curve"]}.pdf'
        LOGGER.info(f'Writing grouped L-curve figure: {pdffile}')
        with PdfPages(pdffile) as pdf:
            # add some info to the PDF
            d=pdf.infodict()
            d['Title']='L-Curve Results'
            d['Author']=pwd.getpwuid(os.getuid()).pw_gecos
            d['Subject']=f'L-Curve results for grouped data from {__code__}'
            d['Keywords']=f'{__code__} WFSS L-curve groups'
            d['Producer']=__code__
            
            for grpid,srcdict in enumerate(groups(sources)):

                self.matrix.build_matrix(data,srcdict)

                res=self.optimizer(self.matrix)
                unc=self.matrix.compute_uncertainty()

                
                for sid,lid in self.matrix.ri.items():
                    segid,regid=self.matrix.sedkeys[sid]
                    lamid=self.matrix.lamids[lid]
                
                    if hasattr(sources[segid],'extpars'):
                        pars=sources[segid].extpars
                    else:
                        pars=self.matrix.defpars

                    # variables for outputing
                    wave=pars.wavelengths()
                    flam=np.full_like(wave,np.nan,dtype=float)
                    func=np.full_like(wave,np.nan,dtype=float)

                    # nota bene:  the methods LSQR and LSMR will formally
                    # return the uncertainties, but are the diagonal
                    # elements of (A'A + ell *I)^-1, which will be
                    # incorrect if ell != 0 (see Ryan Casertano Pirzkal 2018)
                    # Therefore, we also estimate the diagonal elements of
                    # A'A, and invert them.  This is not correct, but a
                    # different error than including the damping.

                    
                    # populate the results
                    flam[lamid]=res.x[lid]
                    func[lamid]=unc[lid]      # from matrix elements
                    #func[lamid]=res.lo[lid]   # from matrix inversion (wrong!)
            
                    # update the outputting sturctures
                    sources[segid].grpid=grpid
                    #sources[segid].spectralregions[regid].sed.reset(wave,flam,
                    sources[segid][regid].sed.reset(wave,flam,func=func)
                                                                    

                # update results for the L-curve data
                kwargs={'nobj':(len(srcdict),'Number of sources')}
                hdulL.append(self.matrix.lcurve.as_HDU(grpid=grpid,**kwargs))
                self.matrix.lcurve.pdfplot(pdf)
            


        # loop over sources for outputting
        for source in sources.values():
        
            hdu=source.as_HDU()            
            if source.is_compound:
                hdul3.append(hdu)
            else:
                hdul1.append(hdu)

            

#        self.matrix.build_matrix(data,sources)
#            
#        res=self.optimizer(self.matrix)
#        
#
#        # loop over elements
#        for sid,lid in self.matrix.ri.items():
#            segid,regid=self.matrix.sedkeys[sid]
#            lamid=self.matrix.lamids[lid]
#            
#            if hasattr(sources[segid],'extpars'):
#                pars=sources[segid].extpars
#            else:
#                pars=self.matrix.defpars
#
#            # variables for outputing
#            wave=pars.wavelengths()
#            flam=np.full_like(wave,np.nan,dtype=float)
#            func=np.full_like(wave,np.nan,dtype=float)
#
#            
#            # populate the results
#            flam[lamid]=res.x[lid]
#            func[lamid]=res.lo[lid]
#            
#            # update the outputting sturctures
#            sources[segid].spectralregions[regid].sed.reset(wave,flam,func=func#)
#
#            # output ascii spectrum
#            #with open(f'multi_{segid}_{regid}.csv','w') as fp:
#            #    print('wavelength,flam,func',file=fp)
#            #    for w,f,df in zip(wave,flam,func):
#            #'        print(f'{w},{f},{df}',file=fp)
#
#            
#
#        # loop over sources for outputting
#        for source in sources.values():
#        
#            hdu=source.as_HDU()            
#            if source.is_compound:
#                hdul3.append(hdu)
#            else:
#                hdul1.append(hdu)
#                
#            
#
#        
#        self.matrix.lcurve.plot(f'{root}_{SUFFIXES["L-curve"]}.pdf')


        # write out the files
        if len(hdul1)>1:
            hdul1[0].header['FILETYPE']='1d spectra'

            x1dfile=f'{root}_{SUFFIXES["1d spectra"]}.fits'
            LOGGER.info(f'Writing 1d extractions: {x1dfile}')
            hdul1.writeto(x1dfile,overwrite=True)
        else:
            LOGGER.warning('No 1d spectra written')
            
        if len(hdul3)>1:
            hdul3[0].header['FILETYPE']='3d spectra'

            x3dfile=f'{root}_{SUFFIXES["3d spectra"]}.fits'
            LOGGER.info(f'Writing 3d extractions: {x3dfile}')
            hdul1.writeto(x3dfile,overwrite=True)
        else:
            LOGGER.warning('No 3d spectra written')


        lcvfile=f'{root}_{SUFFIXES["L-curve"]}.fits'
        LOGGER.info(f'Writing L-curve tabular data {lcvfile}')
        hdulL.writeto(lcvfile,overwrite=True)
