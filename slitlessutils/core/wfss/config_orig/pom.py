import numpy as np
from astropy.io import fits

from ...utilities import headers

def recast(func):
    """
    Custom decorator function to facilitate numpy typing
    
    Parameters
    ----------
    func : callable
        The function to call
    
    Returns
    -------
    inner : callable
        The decorated function
    """

    
    def _inner(self,x,y):
        
        x=np.asarray([x]) if np.isscalar(x) else np.asarray(x)
        y=np.asarray([y]) if np.isscalar(y) else np.asarray(y)

        #x=np.asarray([x]) if np.isscalar(x) else x if isinstance(x,np.ndarray) else np.asarray(x)
        #y=np.asarray([y]) if np.isscalar(y) else y if isinstance(y,np.ndarray) else np.asarray(y)

        p=func(self,x,y)
        if np.size(p,axis=None)==1:
            p=p.item()
        else:
            p=np.squeeze(p)
                
        return p
    return _inner


class POM:
    """
    A base class for implementing a pick-off mirror (POM).  This effectively
    truncates the spatial region around the detector where a potential 
    source may project light (at any wavelength of any order) onto the 
    detector.
    """

    def update_header(self,h):
        """
        Method to update a fits header

        Parameters
        ----------
        h : `astropy.io.fits.Header`
           The fits header to update
        """
        h.set('POM',value=True,comment='Was a POM used?')
        h.set('POMTYPE',value=self.__class__.__name__,comment='POM Type')
        headers.add_stanza(h,'POM Information',after='POM')



class UnityPOM(POM):
    """
    Unity POM.  Returns 1 for all pixels
    """

    def __init__(self):
        pass

    def __call__(self,x,y):
        """ 
        Method to call the POM.

        Parameters
        ---------
        x : float, int, or `np.ndarray`
           The x-coordinates in the dispersed image.

        y : float, int, or `np.ndarray`
           The y-coordinates in the dispersed image.
        
        Returns
        -------
        p : `np.ndarray`
            The POM value.  Will be all `True` for this `UnityPOM`
        
        """

        return np.ones_like(x,dtype=bool)


    def update_header(self,h):
        """
        Method to update a fits header

        Parameters
        ----------
        h : `astropy.io.fits.Header`
           The fits header to update
        """
        super().update_header(h)


class RangePOM(POM):
    """
    Class for implementing a simple rectangular POM

    """

    def __init__(self,xr,yr):
        """
        Initialize a `RangePOM`

        Parameters
        ----------
        xr : tuple, list, `np.ndarray'
            two-element iterable that gives range in x-coordinate

        yr : tuple, list, `np.ndarray'
            two-element iterable that gives range in y-coordinate

        """

        assert len(xr)==2,'XR requires two-element array'
        assert len(yr)==2,'YR requires two-element array'

        self.x0=min(*xr)
        self.x1=max(*xr)
        self.y0=min(*yr)
        self.y1=max(*yr)

        
    @recast
    def __call__(self,x,y):
        """
        Call the POM to test if pixel is in range

        Parameters
        ----------
        x : float, int, or `np.ndarray`
            The x-coordinate in the detector units

        y : float, int, or `np.ndarray`
            The y-coordinate in the detector units

        Returns
        -------
        p : bool or `np.ndarray`
            Returns `True` if the pixel is in range and `False` if not
        
        Notes
        -----
        The input coordinates `x` and `y` must have the same shape, 
        and the returned `p` values will be same shape
        """
        
        return (self.x0<x) & (x<self.x1) & (self.y0<y) & (y<self.y1)

        

    def update_header(self,h):
        """
        Method to update a fits header
        
        Parameters
        ----------
        h : `astropy.io.fits.Header`
           The fits header to update
        
        Returns
        -------
        None
        
        """
        super().update_header(h,used=True)

        h.set('POMX0',value=self.x0,comment='lower-x bound of POM')
        h.set('POMX1',value=self.x1,comment='upper-x bound of POM')
        h.set('POMY0',value=self.y0,comment='lower-y bound of POM')
        h.set('POMY1',value=self.y1,comment='upper-xybound of POM')


class PolygonPOM(RangePOM):
    """
    Class to implement a polygon POM.  

    """

    def __init__(self,px,py):
        """
        Initialize the POM
    
        Parameters
        ----------
        px : list, tuple, or `np.ndarray`
            An interable of the polygon x-coordinates

        py : list, tuple, or `np.ndarray`
            An interable of the polygon y-coordinates
        """
        assert (len(px)>=3),'Must have a x-polygon, which requires >=3 points'
        assert (len(py)>=3),'Must have a y-polygon, which requires >=3 points'

        xr=(np.amin(px),np.amax(px))
        yr=(np.amin(py),np.amax(py))
        RangePOM.__init__(self,xr,yr)

        # gotta have this
        self.px1=np.array(px,dtype=float)
        self.py1=np.array(py,dtype=float)

        # this makes it faster
        self.px2=np.roll(self.px1,-1)
        self.py2=np.roll(self.py1,-1)

    @recast
    def __call__(self,x,y):
        """
        Call the POM to test if pixel is in range

        Parameters
        ----------
        x : float, int, or `np.ndarray`
            The x-coordinate in the detector units

        y : float, int, or `np.ndarray`
            The y-coordinate in the detector units

        Returns
        -------
        p : bool or `np.ndarray`
            Returns `True` if the pixel is in range and `False` if not
        
        Notes
        -----
        The input coordinates `x` and `y` must have the same shape, 
        and the returned `p` values will be same shape
        """


        p=np.zeros_like(x,dtype=bool)

        g=np.where(super().__call__(x,y))[0]
        if len(g)>0:
            
            dx1=self.px1-x[g,np.newaxis]
            dy1=self.py1-y[g,np.newaxis]

            dx2=self.px2-x[g,np.newaxis]
            dy2=self.py2-y[g,np.newaxis]
            

            ang=np.arctan2(dx1*dy2-dy1*dx2,
                           dx1*dx2+dy1*dy2)

            tot=np.abs(np.sum(ang,axis=1))

            p[g]=(tot>np.pi)
            

        
        return p
        
        
    def update_header(self,h):
        """
        Method to update a fits header
        
        Parameters
        ----------
        h : `astropy.io.fits.Header`
           The fits header to update
        
        Returns
        -------
        None
        
        """
        super().update_header(h)
        h.set('POMX',value=pomx,comment='x-coordinates of polygon')
        h.set('POMY',value=pomy,comment='y-coordinates of polygon')
        
        
class ImagePOM(RangePOM):
    """
    Class to implement an image POM (used by JWST)
    """

    def __init__(self,filename,threshold=None,**kwargs):
        """
        Method to initialize the Image POM
        
        Parameters
        ----------
        filename : str
            Full path to a fits file containing the POM info
        
        kwargs : dict, optional
            Optional keywords passed to `astropy.io.fits.getdata()`

        """
        self.filename=filename
        self.data=fits.getdata(self.filename,**kwargs)
        ny,nx=self.data.shape
        RangePOM.__init__(self,(0,nx-1),(0,ny-1))

        self.threshold=threshold
        

    @recast
    def __call__(self,x,y):
        """
        Call the POM to test if pixel is in range

        Parameters
        ----------
        x : float, int, or `np.ndarray`
            The x-coordinate in the detector units

        y : float, int, or `np.ndarray`
            The y-coordinate in the detector units

        Returns
        -------
        p : bool or `np.ndarray`
            Returns `True` if the pixel is in range and `False` if not
        
        Notes
        -----
        The input coordinates `x` and `y` must have the same shape, 
        and the returned `p` values will be same shape
        """
        #raise NotImplementedError('ImagePOM not finished')

        
        #https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.RectBivariateSpline.html#scipy.interpolate.RectBivariateSpline
        
        p=np.zeros_like(x,dtype=float)

        g=np.where(super().__call__(x,y))[0]

        if len(g)>0:
            p[g]=self.data[y[g].astype(int),x[g].astype(int)]

        if self.threshold:
            p>=self.threshold

            
        return p
        


        

    def update_header(self,h):
        """
        Method to update a fits header
        
        Parameters
        ----------
        h : `astropy.io.fits.Header`
           The fits header to update
        
        """
        super().update_header(h)
        hdr.set('POMFILE',value=self.filename,comment='filename of the POM')
        hdr.set('POMTHR',value=self.threshold,comment='threshold for valid')
        hdr.set('POMNX',value=self.x1+1,comment='size in x')
        hdr.set('POMNY',value=self.y1+1,comment='size in y')



def load_pom(**kwargs):
    """
    A helper function to load a POM

    Parameters
    ----------
    kwargs : dict, optional
        key/value pairs to configure which type of POM to return 

    Returns
    -------
    pom : `POM`
        The pick-off mirror result

    Notes
    -----
    if `XRANGE` and `YRANGE` are set, then return `RangePOM`
    elif `POMX` and `POMY` are set, then return `PolygonPOM`
    elif `POMFILE` is set, then return `ImagePOM`
    else return a `UnityPOM`
    """


    if 'XRANGE' in kwargs and 'YRANGE' in kwargs:
        return RangePOM(kwargs['XRANGE'],kwargs['YRANGE'])
    elif 'POMX' in kwargs and 'POMY' in kwargs:
        return PolygonPOM(kwargs['POMX'],kwargs['POMY'])
    elif 'POMFILE' in kwargs:
        return ImagePOM(kwargs['POMFILE'],threshold=kwargs.get('POMTHRESH'))
    else:
        return UnityPOM()
