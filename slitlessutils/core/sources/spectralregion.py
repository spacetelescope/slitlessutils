import numpy as np

from ..photometry import SED


class SpectralRegion:
    """
    class to implement a spectral region.  Here a 'spectral region' is 
    the collection of direct-image pixels that have the same spectrum
    """
    
    def __init__(self,x,y,w,segid,regid,ltv=(0.,0.)):
        """
        Initializer

        Parameters
        ----------
        x : int or `np.ndarray`
           The x-coordinates.  Will be internally recast as `np.ndarray`
        
        y : int or `np.ndarray`
           The y-coordinates.  Will be internally recast as `np.ndarray`

        w : float or `np.ndarray`
           The pixel weights.  Will be internally recast as `np.ndarray`

        segid : int
           The segmentation ID for this region

        regid : int
           The spectral region ID

        ltv : 2-tuple, optional
           The LTV for this positions.  Default is (0,0)

        """

        self.segid=segid
        self.regid=regid
        self.x=np.array(x,dtype=int)
        self.y=np.array(y,dtype=int)
        self.w=np.array(w,dtype=float)
        self.sed=SED()
        self.ltv=ltv

    def __len__(self):
        """
        Method to facilitate the len() operator
        """
        return len(self.x)

    def __iter__(self):
        """
        Method to facilitate an iterator

        Returns
        -------
        x : int
            pixel x-coordinate

        y : int
            pixel y-coordinate

        w : float
            pixel weights
        """
        yield from zip(self.x,self.y,self.w)

    def __getitem__(self,i):
        """
        Method to implement a list-like behavior
        
        Parameters
        ----------
        i : int
           The index to retrieve.  If out-of-bounds, then return `None`
        
        Returns
        -------
        x : int
            pixel x-coordinate

        y : int
            pixel y-coordinate

        w : float
            pixel weights
        """

        try:
            value=(self.x[i],self.y[i],self.w[i])
        except:
            value=None
        return value

    def __str__(self):
        return f'Spectral region: {self.name}'
    
    @property
    def name(self):
        """
        The name of this `SpectralRegion`

        Returns
        -------
        name : 2-tuple
            The name as a tuple of the segid and regid
        """        
        return (self.segid,self.regid)


    def image_coordinates(self,x,y,dtype=int):
        """
        Method to transform coordinates by the LTV keywords
        
        Parameters
        ----------
        x : int, float, or `np.ndarray`
            The x-coordinates

        y : int, float, or `np.ndarray`
            The y-coordinates

        dtype : type or None, optional
            Variable to recast the data types.  If None, then no retyping
            is done.  Default is None

        Returns
        -------
        xx : arbitrary type
            The transformed x-coordinates

        yy : arbitrary type
            The transformed y-coordinates
        
        Notes
        -----
        The input coordinates must have the same shape, their dtype does 
        not matter.  The output coordinates will have that shape.
        """


        xx=x-self.ltv[0]
        yy=y-self.ltv[1]
        if dtype is not None:
            xx=xx.astype(dtype,copy=False)
            yy=yy.astype(dtype,copy=False)

        return xx,yy

    def set_spectral_parameters(self,**kwargs):
        """
        Method to set the spectral parameters to this object as attributes

        Parameters
        ----------
        kwargs : dict, optional
            Dictionary of keywords, can be 'wave0','wave1','dwave','scale'

        """
        
        parnames=('wave0','wave1','dwave','scale')
        for k,v in kwargs.items():
            if k in parnames:
                setattr(self,k,v)

    def get_spectral_parameters(self):
        """
        Method to get the spectral parameters from this object's attributes

        Returns
        -------
        pars : tuple
            A tuple of the extraction parameters and the seg/reg IDs
        
        """
        
        pars=[self.segid,self.regid]
        for parname in ('wave0','wave1','dwave','scale'):
            if hasattr(self,parname):
                pars.append(float(getattr(self,parname)))
            else:
                pars.append(np.nan)
        pars=tuple(pars)
        
        #pars={}
        #for parname in ('wave0','wave1','dwave','scale'):
        #    if hasattr(self,parname):
        #        pars[parname]=getattr(self,parname)
                
        return pars
    
    def pixels(self,dtype=int,weights=False,applyltv=False,**kwargs):
        """
        A generator to loop over all direct image pixels

        Parameters
        ----------
        applyltv : bool, optional
            Flag to apply the LTV before returning.  Internally to this 
            source, the pixels are stored in lowest integers, and not on
            the original pixel grid.  Default is False
        
        dtype : type, optional
            The dtype to return the coordinates.  Default is int.

        Returns
        -------
        x : arb type
            The x-coordinates
        
        y : arb type
            The y-coordinates

        """


        if applyltv:
            if weights:
                for x,y,w in zip(self.x,self.y,self.w):
                    xx,yy=self.image_coordinates(x,y,**kwargs)
                    yield xx,yy,w
            else:
                for x,y in zip(self.x,self.y):
                    yield self.image_coordinates(x,y,**kwargs)
        else:
            if weights:
                yield from zip(self.x,self.y,self.w)
            else:
                yield from zip(self.x,self.y)
