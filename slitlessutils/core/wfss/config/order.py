import numpy as np
from astropy.io import fits

from ...photometry import Band
#from .sensitivity import Sensitivity
from ....logger import LOGGER


#
# Written for slitlessutils by R. Ryan at STScI
#
# v1.0 Jan 2023
#
#
#




class Order:
    """
    Class to contain the trace, dispersion, and sensitivity for a spectral
    order.
    """

    def __init__(self,order):
        """
        Initialize the order class.

        Parameters
        ----------
        order : str
            The name of the spectral order.  Should be a string, but is 
            retyped to string regardless.

        Returns
        -------
        None
        """


        self.order=str(order)
        self.sensitivity=lambda x: np.full_like(x,np.nan)
        self.dispx=ParametricPolynomial('dispx')
        self.dispy=ParametricPolynomial('dispy')
        self.displ=ParametricPolynomial('displ')

    def __str__(self):
        s=f'Configuration for spectral order: {self.order}\n'
        return s
        
    def load_sensitivity(self,sensfile,**kwargs):
        """
        Method to load the sensitivity curve from a filename.

        Parameters
        ----------
        sensfile : str
            Full path to the sensitivity file.

        kwargs : dict, optional 
            keywords that get passed to the `Sensitivity` curve object.
        """

        self.sensitivity=Sensitivity(sensfile,**kwargs)

    def dispersion(self,x0,y0,wavelength=None):
        """
        Method to compute the dispersion in A/pix at some undispersed
        position and wavelength.  This is given by the derivative of the 
        wavelength solution as a function of position along the spectral 
        trace.

        .. math::
           \frac{d\lambda}{dr} = \frac{d\lambda/dt}{\sqrt{(dx/dt)^2 + (dy/dt)^2}}
        where :math:`t` is parametric value given by:
        .. math::
           t = DISPL^{-1}(\lambda)
        
        and :math:`DISPL` comes from the grismconf file.

        Parameters
        ----------
        x0 : int or float
           The undispersed x-coordinate

        y0 : int or float
           The undispersed y-coordinate

        wavelength : int, float, or None, optional
           The wavelength (in A) to compute the dispersion.  If set to 
           `None`, then the bandpass-averaged wavelength is used.
        
        Returns
        -------
        dldr : float
           The instaneous dispersion in A/pix.

        """


        # default wavelength to the average value over the sensitivity
        if wavelength is None:
            wavelength=self.sens.wave_ave

        # compute the dispersion using:
        # t = h^-1(lambda)
        # x(t) = a+ b*t + c*t^2 + ...
        # x'(t) = b+ 2*c*t + ...
        # y'(t) = r+ 2*s*t + ...
        # l'(t) = u+ 2*v*t + ...
        # r'(t) = sqrt(dxdt**2 + dydt**2)
        # dl/dr = l'(t)/r'(t)

        # invert the DISPL function
        t=self.displ.invert(xyg,wavelength)

        # evaluate the derivatives (w.r.t. to t) of the DISPX,
        # DISPY, and DISPL functions, which came from grism conf
        dxdt=self.dispx.deriv(xyg,t)
        dydt=self.dispy.deriv(xyg,t)
        dldt=self.displ.deriv(xyg,t)

        # return the dispersion
        dldr=dldt/np.sqrt(dxdt*dxdt+dydt*dydt)
        return dldr
    

    def deltas(self,x0,y0,wav):
        """
        Method to compute the offsets  with respect to the undispersed 
        position.  NOTE: the final WFSS position would be given by 
        adding back the undispersed positions.

        Parameters
        ----------
        x0 : float, `np.ndarray`, or int
           The undispersed x-position.
        
        y0 : float, `np.ndarray`, or int
           The undispersed y-position.
        
        wav : `np.ndarray`
           The wavelength (in A).

        Returns
        -------
        dx : float or `np.ndarray`
           The x-coordinate along the trace with respect to the undispersed 
           position.

        dy : float or `np.ndarray`
           The y-coordinate along the trace with respect to the undispersed 
           position.

        Notes
        -----
        The undispersed positions (`x0` and `y0`) must be of the same 
        shape, hence either both scalars or `np.ndarray`s with the same 
        shape.  The dtype of the variables does not matter.  If these 
        variables are arrays, then the output will be a two-dimensional 
        array of shape (len(wav),len(x0)).  If they are scalars, then
        the output will be a one-dimensional array with shape (len(wav),).

        """
        
        if np.isscalar(x0):
            xy0=(x0,y0)
            t=self.displ.invert(xy0,wav)
            dx=self.dispx.evaluate(xy0,t)
            dy=self.dispy.evaluate(xy0,t)
        else:
            nx=len(x0)
            nw=len(wav)
            # fill the data
            dx=np.empty((nw,nx),dtype=float)
            dy=np.empty((nw,nx),dtype=float)
            for i,xy0 in enumerate(zip(x0,y0)):
                t=self.displ.invert(xy0,wav)
                dx[:,i]=self.dispx.evaluate(xy0,t)
                dy[:,i]=self.dispy.evaluate(xy0,t)

        return dx,dy

    
    #def disperse(self,x0,y0,wav):
    #    dx,dy=self.deltas(x0,y0,wav)
    #    return dx+x0,dy+y0
         

    @property
    def name(self):
        """
        The name of the order
        """
        
        return str(self.order)






class Sensitivity(Band):
    """
    Class that holds the WFSS sensitivity curve. 
    
    inherits from `su.core.photometry.Band`
    """
    
    def __init__(self,sensfile,senslimit=1e10):
        """
        Initialize the `Sensitivity` object.
        
        Parameters
        ----------
        sensfile : str
            The full path to the sensitivity file, which should be a fits
            table.

        senslimit : float
            The limit, below which, the data are ignored for the purposes
            of computing the average values.  Default is 1e10.

        """
        self.sensfile=sensfile
        self.senslimit=senslimit

        # read the file
        data,header=fits.getdata(self.sensfile,1,header=True)

        g=np.where(data['SENSITIVITY'] > self.senslimit)
        
        
        Band.__init__(self,data['WAVELENGTH'],data['SENSITIVITY'],
            where=g,unit=header.get('TUNIT1',''))


    def __str__(self):
        return f'Sensitivity curve: {self.sensfile}'

  
class ParametricPolynomial(list):
    """
    Class to implement the parametric polynomial.

    .. math::
       f(t;x_0,y_0) = a(x_0,y_0) + b(x_0,y_0) t + c(x_0,y_0) t^2 + ...
    
    where each coefficient is given as a `SpatialPolynomial`.

    inherits from list
    """
    
    def __init__(self,name):
        """
        Initialize the object.

        Parameters
        ----------
        name : str
            The name of the polynomial.
        """
        

        self.name=name
        self.order=-1
        
        # provide a fallback
        self.invert = lambda xy,f: None

    def append(self,coefs):
        """
        Method to include another parametric order (overrides the 
        list.append method).
        
        Parameters
        ----------
        coefs : list, tuple, or `np.ndarray`
            The coefficients to pass to the `SpatialPolynomial` object.
            The length of this iterable must be a triangular number.
        
        """

        poly=SpatialPolynomial(coefs)
        if poly:
            super().append(poly)
            self.order+=1
        
            if self.order==1:
                self.invert=self._first
            elif self.order==2:
                self.invert=self._second
            else:
                self.invert=self._nth

    def _first(self,xy,f):
        """
        Helper method to analytically invert a 1st order polynomial
        
        Parameters
        ----------
        xy : tuple, list, or `np.ndarray`
           The (x,y) pair.  
        
        f : float
           The value of the polynomial
        
        Returns
        -------
        t : float
           The parameter that gives the value `f`.
        """

        coefs=self.coefs(xy)
        t=(f-coefs[0])/coefs[1]
        return t

    def _second(self,xy,f):
        """
        Helper method to analytically invert a 2nd order polynomial
        
        Parameters
        ----------
        xy : tuple, list, or `np.ndarray`
           The (x,y) pair.  
        
        f : float
           The value of the polynomial
        
        Returns
        -------
        t : float
           The parameter that gives the value `f`.
        """

        coefs=self.coefs(xy)


        disc=coefs[1]*coefs[1]-4.0*(coefs[0]-f)*coefs[2]

        sqrt=np.sqrt(disc)

        
        nump=(-coefs[1]+sqrt)/(2*coefs[2])
        #numm=(-coefs[1]-sqrt)/(2*coefs[2])

        return nump


    def _nth(self,xy,f):
        ''' numerically invert an arbitrary order polynomial '''
        raise NotImplementedError('arbitrary order not implemented.')

    
    def coefs(self,xy):
        """ 
        Method to evaluate all the `SpatialPolynomial` coefficients
        for a given position.

        Parameters
        ----------
        xy : list, tuple, `np.ndarray`
            the (x,y) pair for which the list of `SpatialPolynomial`s 
            are evaluated.
        
        Returns
        -------
        coefs : list
            The coefficients for each parameteric order 

        """

        coefs=[poly.evaluate(xy) for poly in self]

        #poly=np.poly1d(coefs[::-1])
        
        return coefs

    def evaluate(self,xy,t):
        """
        Method to evaluate the polynomial at some position and parameter
        
        Parameters
        ----------
        xy : list, tuple, `np.ndarray`
            The spatial positions to pass to the `SpatialPolynomial`s

        t : float, int, or `np.ndarray`
            The parameter to evaluate this `ParametricPolynomial`

        Returns
        -------
        f : float, or `np.ndarray`
            The value of the polynomial
        """
        
        f=sum(p.evaluate(xy)*t**i for i,p in enumerate(self))
        return f
        
    def deriv(self,xy,t):
        """
        Method to evaluate the derivative of the polynomial at some 
        position and parameter
        
        Parameters
        ----------
        xy : list, tuple, `np.ndarray`
            The spatial positions to pass to the `SpatialPolynomial`s

        t : float, int, or `np.ndarray`
            The parameter to evaluate this `ParametricPolynomial`

        Returns
        -------
        dfdt : float, or `np.ndarray`
            The value of the derivative of the polynomial
        """
        
        dfdt=sum(p.evaluate(xy)*i*t**(i-1) for i,p in enumerate(self[1:],start=1))
        return dfdt 
  





class SpatialPolynomial(dict):
    """
    Class to implement 2d spatial polynomial of the form:
    
    .. math::
    
       p(x,y) = a + bx + cy + dx^2 + exy + fy^2 + ...


    inherts from dict.  The key/value pairs are the spatial exponents
    (stored via 'Cantor pairing') and coefficients, respectively.

    """
    
    # used for printing
    EXP={0:'\u2070',1:'\u00b9',2:'\u00b2',3:'\u00b3',4:'\u2074',
         5:'\u2075',6:'\u2076',7:'\u2077',8:'\u2078',9:'\u2079'}
    

    def __init__(self,values):
        """
        Method to initialize the object.

        Parameters
        ----------
        values : float, int, or `np.ndarray`
            The coefficients for this `SpatialPolynomial`.  Due to the 
            definition of the spatial polynomial, the number of 
            elements of `values` must be a triangular number.  See:
            
            https://en.wikipedia.org/wiki/Triangular_number
        
        """

        
        self.order=None
        if np.isscalar(values):
            self.order=0
            self[(0,0)]=values
        else:
            # the coefs are stored as cantor pairs.
            # https://en.wikipedia.org/wiki/Pairing_function#Inverting_the_Cantor_pairing_function]

            n=(np.sqrt(1+8*len(values))-1)/2
            if n.is_integer():
                n=int(n)
                self.order=n-1

                # old way of decoding cantor pairing
                i=0
                for j in range(n):
                    for k in range(j+1):
                        self[(j-k,k)]=values[i]
                        i+=1
            else:
                msg="Input must be an array whose length is a triangular number"
                LOGGER.error(msg)
                raise RuntimeError(msg)
                

            
    def __str__(self):
        s=f"Spatial polynomial of order {self.order}:\n"

        p=[f'{cij} x{self.EXP[i]} y{self.EXP[j]}' for (i,j),cij in self.items()]
        s+='p(x,y)='+'+'.join(p)
        return s

    def evaluate(self,xy):
        """
        Method to evaluate this spatial polynomial at a position

        Parameters
        ----------
        xy : list, tuple, `np.ndarray`
            An iterable with 2 elements, where the two elements are for 
            the `x` and `y`.  
    
        Returns
        -------
        p : float, `np.ndarray`
            The value of the polynomial
        """
        p=sum(coef*xy[0]**i*xy[1]**j for (i,j),coef in self.items())
        return p

if __name__=='__main__':
    x=Order(2)
    print(x)
