import numpy as np
from scipy import ndimage
from skimage import measure
import ds9regions as ds9reg

from ...utilities import indices
from ..module import Module
from ....logger import LOGGER


class Region(Module):
    """
    Module to create ds9-style regions files
    
    inherits from `su.core.modules.Module`

    Parameters
    ----------
    orders : list or None, optional
       The spectral orders to make regions.  If `None`, then will do all 
       orders present in the configuration file.  Default is None
       
    nthin : int, optional
       The thinning frequency of points to put in the ds9 region, default is 4

    close_size : int, optional
       The size of the morphological closing operator, default is 9
       
    kwargs : dict, optional
       Optional keywords sent to parent class, see `su.core.modules.Module`
    
    """


    DESCRIPTION = 'Regions'
    
    COLORS={'0':'white','+1':'#1f77b4','+2':'#ff7f0e','+3':'#2ca02c',
            '+4':'#d62728','+5':'#9467bd','+6':'#8c564b','+7':'#e377c2',
            '+8':'#7f7f7f','+9':'#bcbd22','+10':'#17becf','-1':'#aec7e8',
            '-2':'#ffbb78','-3':'#98df8a','-4':'#ff9896','-5':'#c5b0d5',
            '-6':'#c49c94','-7':'#f7b6d2','-8':'#c7c7c7','-9':'#dbdb8d',
            '-10':'#9edae5'}



    def __init__(self,orders=None,nthin=4,close_size=9,**kwargs):
        Module.__init__(self,self.region,**kwargs)
        self.orders=orders
        self.close_size=close_size
        self.nthin=nthin

        
    def region(self,data,sources):
        """
        Function to create the regions
        
        Parameters
        ----------
        data : `WFSSCollection`
           The WFSS data to tabulate

        sources : `SourceCollection`
           The sources to tabulate
       
        Returns
        -------
        result : tuple
           A tuple of the region files created

        Notes
        -----
        There will be as many files as there are detectors in one of the 
        elements of the data.  They will have the following filenames:
        
        >>> f'{data.dataset}_{detdata.name}.reg'

        
        """


        # an internal variable that doesn't really affect things much
        pad=7
        
        # a morphologoical closing operator
        struct=np.ones((self.close_size,self.close_size),dtype=float)

        # a bunch of output files        
        outfiles=[]
        
        # get wavelengths
        #wav=data.grating.wavelengths(nsub=1./2.)
        wav=data.disperser.wavelengths(nsub=2)

        kwargs={'width':4,'move':False,'rotate':False,'fixed':True,'edit':False}
        for detname,detdata in data.items():
            naxis=detdata.naxis      # this just saves typing later
            orders=self.orders if self.orders else detdata.orders


            regfile=f'{data.dataset}_{detdata.name}.reg'
            outfiles.append(regfile)
            with ds9reg.DS9Regions(regfile,**kwargs) as rf:
                            
                for segid,source in sources.items():
                    x,y=[],[]
                    for args in source.pixels():
                        x.append(args[0])
                        y.append(args[1])
                    x=np.asarray(x)
                    y=np.asarray(y)

                    xx,yy=detdata.xy2xy(x,y,source.wcs,forward=False)
                    

                    for order in orders:
                        xg,yg=detdata.config.config.disperse(xx,yy,order,wav)
                        
                        
                                                
                        xg=xg.astype(int).flatten()
                        yg=yg.astype(int).flatten()



                        x0=np.amin(xg)-pad
                        x1=np.amax(xg)+pad
                        y0=np.amin(yg)-pad
                        y1=np.amax(yg)+pad
                        
                        # image size
                        dim=(y1-y0+1,x1-x0+1)

                        
                        yxg=np.ravel_multi_index((yg-y0,xg-x0),dim,order='F')
                        yxg=indices.uniq(yxg)
                        yg,xg=np.unravel_index(yxg,dim,order='F')
                        xg+=x0
                        yg+=y0

                        #idx=np.ravel_multi_index((xg,yg),dim,order='F')
                        #idx=indices.uniq(idx)
                        #xg,yg=np.unravel_index(idx,dim,order='F')


                        # fill it one                        
                        msk=np.zeros(dim,dtype=int)
                        #msk[yg,xg]=1
                        msk[yg-y0+1,xg-x0+1]=1
                        msk=ndimage.binary_closing(msk, structure=struct)
                        
                        # get the contour
                        contours=measure.find_contours(msk.astype(int),0.1,
                                                       fully_connected='high')
                        # make the contours and thin
                        xc=contours[0][::self.nthin,1]+x0
                        yc=contours[0][::self.nthin,0]+y0

                        
                        
                        # create the region
                        color=self.COLORS.get(order,'green')
                        reg=ds9reg.Polygon(xc,yc,color=color)

                        # write to the file
                        rf.write_region(reg)
        return tuple(outfiles)
