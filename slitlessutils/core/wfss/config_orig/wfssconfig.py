import numpy as np
import os

from .flatfield import load_flatfield
from .order import Order
from .pom import load_pom 


#
# Written for slitlessutils by R. Ryan at STScI
#
# v1.0 Jan 2023
#
#



class WFSSConfig(dict):
    """
    Primary class to configure a WFSS.  Contains information on the 
    trace/dispersion, flatfield, pick-off-mirror (POM), and sensitivity.
    
    Inherits from dict, where each keyword,value pair is a spectral order
    name and order `slitlessutils.core.wfss.config.order.py`, respectively.    

    """
    


    COMMENTS=('#','%','!',';','$')    # comments in the file

    
    def __init__(self,conffile,fwcpos=None,band=None):
        """
        
        Parameters
        ----------
        conffile : str
            Filename for the configuration file.  This should be in 
            the `grismconf` format.

        fwcpos : float or None, optional
            The 'filter wheel position', which is an angle (in degree) 
            of the filter wheel.  Used only for JWST/NIRISS.  A value of 
            'None' is interpreted as 0 deg.  Default value is None.

        band : str or None, optional
            The name of the direct image filter, which may dictate any 
            wedge offsets.  Used only for WFC3/IR.  A value of 'None' is
            interpreted as blank.  Default value is None.         

        """
        
        self.confpath=os.path.dirname(conffile)
        
        self.conffile=conffile

        self.ffname=None
        self.set_rotation(0.)
        self.set_shift(0.,0.)
        
        pomdata={}
        self.data=[]

        
        with open(self.conffile,'r') as fp:
            for line in fp:
                line=line.strip()
                self.data.append(line)
                tokens=line.split()
                if tokens and tokens[0][0] not in self.COMMENTS:
                    key=tokens.pop(0)

                    value=list(map(self.retype,tokens))

                    nval=len(value)
                    if nval==0:
                        order=key.split('_')[1]
                        self[order]=Order(order)
                    else:
                        if nval==1:
                            value=value[0]
                        else:
                            value=np.array(value)

                        # chekc the key for being a valid order info
                        if key == 'FWCPOS_REF' and fwcpos is not None:
                            self.set_rotation(np.radians(value-fwcpos))
                        elif key==f'WEDGE_{band}'.upper() and nval==2:
                            self.set_shift(value[0],value[1])
                        elif key=='FFNAME':
                            self.ffname=value
                        elif key=='XRANGE':
                            pomdata[key]=value
                        elif key=='YRANGE':
                            pomdata[key]=value
                        elif key=="POMX":
                            pomdata[key]=value
                        elif key=='POMY':
                            pomdata[key]=value
                        elif key=='POMFILE':
                            pomdata[key]=value                            
                        else:
                            self._set_order_data(key,value)


        # precisely which type of pom do we have?
        self.pom=load_pom(**pomdata)

    def set_rotation(self,theta):
        """
        Method to set a rotation in the spectral trace.
        
        Parameters
        ----------
        theta : float
            Rotation angle in radians

        Returns
        -------
        None
        """

        
        self.theta=theta
        self.cos=np.cos(theta)
        self.sin=np.sin(theta)

    def add_rotation(self,dtheta):
        """
        Method to update the current spectral dispersion  rotation angle
        
        Parameters
        ----------
        theta : float
            Rotation angle in radians

        Returns
        -------
        None
        """
        self.set_rotation(self.theta+dtheta)        
        
    def set_shift(self,dx,dy):
        """
        Method to set a shift in the spectral trace.
        
        Parameters
        ----------
        dx : float
            Shift in x in pixels.

        dy : float
            Shift in y in pixels.

        Returns
        -------
        None
        """        
        self.dx=dx
        self.dy=dy

    def add_shift(self,dx,dy):
        """
        Method to update the shift in the spectral trace.
        
        Parameters
        ----------
        dx : float
            Shift in x in pixels.

        dy : float
            Shift in y in pixels.

        Returns
        -------
        None
        """        

        self.dx+=dx
        self.dy+=dy
    
        
    def load_flatfield(self,unity=False):
        """
        Method to read a flat field from the configuration file.

        Parameters
        ----------
        unity : boolean, optional
            Flag to return a unity flat.  Default is False

        Returns
        -------
        flatfield : `np.ndarray`
            The two-dimensional flat field
        """
        
        if unity or (self.ffname is None):
            ff=load_flatfield()
        else:
            ff=load_flatfield(self.ffname)
        return ff


    def dispersion(self,x0,y0,order,wavelength=None):
        """
        Method to determine the dispersion (ie. the instaneous derivative 
        of the wavelength as a function of position along the trace) at 
        a given position and spectral order.

        See also `slitlessutils.core.wfss.config.order.py`.

        Parameters
        ----------
        x0 : float
            The x-position for which to compute the dispersion.

        y0 : float
            The y-position for which to compute the dispersion.

        order : str
            The name of the spectral order
        
        wavelength : float or None, optional
            The wavelength to compute the dispersion.  If set as 'None',
            then interpreted as the 'PHOTPLAM' for the spectroscopic 
            sensitivity.  Default is None.
        
        Returns
        -------
        dispersion : float
            The dispersion in A/pix.

        """
        return self[order].dispersion(x0,y0,wavelength=wavelength)


    def disperse(self,x0,y0,order,wav):
        """
        Method to convert an undispersed position into a WFSS image 
        position using the calibrated trace and dispersion polynomials.

        Parameters
        ----------
        x0 : float or `np.ndarray`
            The undispersed x-position (in pixels).

        y0 : float or `np.ndarray`
            The undispersed y-position (in pixels).
        
        order : str
            The name of the spectral order

        wav : float or `np.ndarray`
            The wavelengths (in A) to disperse.  

        Returns
        -------
        x : float or `np.ndarray`
            The WFSS image pixel coordinate.

        y : float or `np.ndarray`
            The WFSS image pixel coordinate.
        
        Notes
        -----
        The coordinates `x0` and `y0` must be of the same datatype. For 
        example, if one is a `np.ndarray`, then they both must be 
        `np.ndarray` (and of the same shape, but dtype is not relevant).
        """

        dx,dy=self[order].deltas(x0,y0,wav)

        x=self.cos*dx-self.sin*dy+self.dx+x0
        y=self.sin*dx+self.cos*dy+self.dy+y0

        return x,y

    def sensitivity(self,order,wav,**kwargs):
        """
        Method to compute the sensitivity for a given order at a 
        selection of wavelengths.

        Parameters
        ----------
        order : str
            The spectral order name

        wav : float or `np.ndarray`
            The wavelengths for which the sensitivity is computed

        kwargs : dict, optional
            Keywords passed to `slitlessutils.core.wfss.config.order.py`,
            which largely govern the interpolation.
        
        Returns
        -------
        sensitivity : float or `np.ndarray`
            The sensitivity (in units of the sensitivity file).  Will be of 
            same datatype and/or shape as `wav`.

        """

        
        return self[order].sensitivity(wav,**kwargs)

    def _set_order_data(self,key,value):
        """
        Helper function in assigning a keyword/value pair to the self.

        This likely should not be used from the outside.

        Parameters
        ----------
        key : str
            The keyword from the grismconf file.
        
        value : float, int, `np.ndarray`, or str
            The value from the grismconf file.
        
        Returns
        -------
        None
        """
        

        tokens=key.split('_')
        n=len(tokens)
        if n==3:
            a=tokens[0] in ('DISPX','DISPY','DISPL')
            b=tokens[1] in self
            c=tokens[2].isdigit()
            if a and b and c:
                getattr(self[tokens[1]],tokens[0].lower()).append(value)
        elif n==2:
            if tokens[0]=='SENSITIVITY' and tokens[1] in self:
                sensfile=os.path.join(self.confpath,value)
                self[tokens[1]].load_sensitivity(sensfile)
            
    @staticmethod
    def retype(datum):
        """
        Static Method to recast a scalar string into either int or 
        float as possible.

        If datum can be an int, then return it. If it can't then try 
        typing it as a float.  If still can't be typed as a float, then
        return the input data.
        
        
        Parameters
        ----------
        datum : str
            Input datum to retype.  Must be a single value.

        Returns
        -------
        datum : int or float or str
            Output retyped datum
        """
        
        try:
            datum=int(datum)
        except ValueError:
            try:
                datum=float(datum)
            except ValueError:
                pass
        return datum

                
if __name__=='__main__':
    x=GrismConfig('/Users/rryan/slitlessutils_config/instruments/WFC3IR/g102.conf')
