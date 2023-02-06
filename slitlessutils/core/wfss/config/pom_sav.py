import numpy as np
from astropy.io import fits

from ...utilities import headers



class POM:
    """
    A base class for implementing a pick-off mirror (POM).  This effectively
    truncates the spatial region around the detector where a potential 
    source may project light (at any wavelength of any order) onto the 
    detector.
    """

    def __init__(self,x0,x1,y0,y1):
        """
        The class initializer
        
        Parameters
        ----------
        x0 : float or int
           The minimum x-coordinate.

        x1 : float or int
           The maximum x-coordinate

        y0 : float or int
           The minimum y-coordinate.

        y1 : float or int
           The maximum y-coordinate
        """

        self.x0=x0
        self.y0=y0
        self.x1=x1
        self.y1=y1


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

        
        #return self.x0<x<self.x1 and self.y0<y<self.y1
        return (self.x0<x) & (x<self.x1) & (self.y0<y) & (y<self.y1)
        
    def update_header(self,h,used=False):
        """
        Method to update a fits header
        
        Parameters
        ----------
        h : `astropy.io.fits.Header`
           Fits header to update

        used : bool, optional
           Flag that a POM was used. Default is False

        """
        
        h['POM']=(False,'Was POM used?')
        h.set('POMTYPE',value=type(self).__name__,comment='Type of POM',
              after='POM')
        h['POMX0']=(self.x0,'lower-bound of x')
        h['POMY0']=(self.y0,'lower-bound of y')
        h['POMX1']=(self.x1,'upper-bound of x')
        h['POMY1']=(self.y1,'upper-bound of y')
        
        headers.add_stanza(h,'POM Settings',before='POM')
        


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
        super().update_header(h,used=True)

        
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

        
        POM.__init__(self,xr[0],xr[1],yr[0],yr[1])


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

        
class PolygonPOM(POM):
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
        raise NotImplementedError('The __call__ method is not ready')
        
        x0,x1=np.amin(px),np.amax(px)
        y0,y1=np.amin(py),np.amax(py)
        POM.__init__(self,x0,x1,y0,y1)

        self.px1=np.array(px,dtype=float)
        self.py1=np.array(py,dtype=float)
        self.px2=np.roll(self.px1,-1)
        self.py2=np.roll(self.py1,-1)

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
        
        if super().__call__(x,y):

            x1=self.px1-x
            y1=self.py1-y
            x2=self.px2-x
            y2=self.py2-y

            dp=x1*x2+y1*y2
            cp=x1*y2-y1*x2
            ang=np.arctan2(cp,dp)

            tot=np.sum(ang)

            return np.abs(tot)>np.pi
        else:
            return np.zeros_like(x,dtype=bool)

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
        h['POMX']=(','.join(str(x) for x in self.px1),'x-coordinates of polygon')
        h['POMY']=(','.join(str(y) for y in self.py1),'y-coordinates of polygon')
        
class ImagePOM(POM):
    """
    Class to implement an image POM (used by JWST)
    """

    def __init__(self,filename,**kwargs):
        """
        Method to initialize the Image POM
        
        Parameters
        ----------
        filename : str
            Full path to a fits file containing the POM info
        
        kwargs : dict, optional
            Optional keywords passed to `astropy.io.fits.getdata()`

        """
        raise NotImplementedError("Image POM are not ready")
        self.filename=filename
        self.data=fits.getdata(self.filename,**kwargs)
        ny,nx=self.data.shape
        POM.__init__(self,0,nx,0,ny)

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
        raise NotImplementedError("__call__ not ready")

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
        h['POMFILE']=(self.filename,'filename of the POM')
        
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
        return ImagePOM(kwargs['POMFILE'])
    else:
        return UnityPOM()

        


        
if __name__=='__main__':
    x=POM.from_ranges([1.,2],[1,2])
    print(x(1.5,1.5))

    y=UnityPOM()
    print(y(2,3))
