import os
import numpy as np

from .asciiconfig import AsciiConfig
from .parametricpolynomial import ParametricPolynomial
from .sensitivity import Sensitivity
from ....logger import LOGGER


class Order:
    ''' class to implement the orders as described by Pirzkal '''

    @classmethod
    def from_conffile(cls,conffile,order):
        conf=AsciiConfig(conffile)
        obj=cls().from_dict(conf,order)
        obj.conffile=conffile
        return obj
    
    @classmethod
    def from_dict(cls,config,order):
        #obj=super(Order,cls).__new__(cls)
        obj=cls()
        obj.conffile=None
        
        obj.order=order

        # three main functions in Pirzkal
        obj.dispx=ParametricPolynomial('dispx')
        obj.dispy=ParametricPolynomial('dispy')
        obj.displ=ParametricPolynomial('displ')

        # support the wedge offsets
        obj.wedges={}

        # extract the data from the config dict
        for key,val in config.items():
            if val is not None:
                nval=1 if np.isscalar(val) else len(val)
            
                if key.startswith(f'XRANGE_{obj.order}') and nval==2:
                    obj.xr=np.array(val,dtype=np.float32)

                elif key.startswith(f'YRANGE_{obj.order}') and nval==2:
                    obj.yr=np.array(val,dtype=np.float32)

                elif key.startswith(f'DISPX_{obj.order}'):
                    obj.dispx.append(val)
                    
                elif key.startswith(f'DISPY_{obj.order}'):
                    obj.dispy.append(val)

                elif key.startswith(f'DISPL_{obj.order}'):
                    obj.displ.append(val)
                
                #elif key.startswith('NAXIS') and nval==2:
                #    obj.naxis=np.array(val,dtype=np.uint16)
                    
                elif key.startswith(f'SENSITIVITY_{obj.order}'):
                    sensfile=os.path.join(config.confpath,val)
                    obj.sensitivity=Sensitivity(sensfile)
                    
                elif key.startswith('WEDGE') and nval==2:
                    filt=key.split('_')[1]
                    obj.wedges[filt]=np.array(val,dtype=np.float32)
                else:
                    pass


        # load the polygon clipping
        #obj.polyclip=polyclip.Polyclip(obj.naxis)

        return obj
    def __bool__(self):
        return bool(self.dispx) and \
            bool(self.dispy) and \
            bool(self.displ)


    @property
    def name(self):
        return self.order

    def __str__(self):
        return f'WFSS order configuration for {self.order}'

    

    def dispersion(self,xyg,wavelength=None):
        ''' compute the dispersion in A/pix at some position and wavelength '''

        # make the default wavelength to the average value over the
        # sensitivity curve
        if wavelength is None:
            wavelength=self.sensitivity.wave_ave
        
        # t = h^-1(lambda)
        # x(t) = a+ b*t + c*t^2 + ...
        # x'(t) = b+ 2*c*t + ...
        # y'(t) = r+ 2*s*t + ...
        # l'(t) = u+ 2*v*t + ...
        # r'(t) = sqrt(dxdt**2 + dydt**2)
        # dl/dr = l'(t)/r'(t)

        t=self.displ.invert(xyg,wavelength)
        dxdt=self.dispx.deriv(xyg,t)
        dydt=self.dispy.deriv(xyg,t)
        dldt=self.displ.deriv(xyg,t)


        #raise RuntimeError("THIS IS NOT CORRECT IN GENERAL")
        #dxdt=self.dispx.coefs(xyg,order=1)
        #dydt=self.dispy.coefs(xyg,order=1)
        #dldt=self.displ.coefs(xyg,order=1)

        dldr=dldt/np.sqrt(dxdt*dxdt+dydt*dydt)
        return dldr

    def wavelengths(self,xc,yc,subsample=1):
        ''' compute the wavelengths at some position '''
        disp=self.dispersion((xc,yc))/subsample
        wave=np.arange(self.sensitivity.wmin,self.sensitivity.wmax,disp,
                       dtype=float)
        return wave

    def disperse(self,xd,yd,wav,band=None):
        ''' compute dispersed positions given a position and wavelengths '''

        # get dimensionalities
        nw=len(wav)
        if np.isscalar(xd):
            xyd=(xd,yd)
            t=self.displ.invert(xyd,wav)
            xg=self.dispx.evaluate(xyd,t)+xyd[0]
            yg=self.dispy.evaluate(xyd,t)+xyd[1]
        else:
            nx=len(xd)

            # fill the data
            xg=np.empty((nw,nx),dtype=float)
            yg=np.empty((nw,nx),dtype=float)
            for i,xyd in enumerate(zip(xd,yd)):
                t=self.displ.invert(xyd,wav)
                xg[:,i]=self.dispx.evaluate(xyd,t)+xyd[0]
                yg[:,i]=self.dispy.evaluate(xyd,t)+xyd[1]

        # apply the wedge offset
        if band is not None:
            if band in self.wedges:
                # NOR HAS THIS A NEGATIVE SIGN
                xg-=self.wedges[band][0]
                yg-=self.wedges[band][1]
            else:
                LOGGER.warning(f'{band} not found in wedge table')

        return xg,yg




    def drizzle(self,xd,yd,wav,ignore='average',band=None,pixfrac=1.):
        ''' compute the fracional pixel area '''

        raise RuntimeError("This gets moved to tabulate")
        
        # apply the ignoring procedure
        if hasattr(self,'xr') and hasattr(self,'yr'):
            ignore=ignore.lower()
            if ignore=='average':
                # if average of pixel is in bounding box
                xave=np.average(xd)
                yave=np.average(yd)
                if (xave<self.xr[0]) or (xave>self.xr[1]) or \
                   (yave<self.yr[0]) or (yave>self.yr[1]):
                    return [],[],[],[]
            elif ignore=='minmax':
                # test min/max in range
                x0,x1=np.amin(xd),np.amax(xd)
                y0,y1=np.amin(yd),np.amax(yd)
                if (x1<self.xr[0]) or (x0>self.xr[1]) or \
                   (y1<self.yr[0]) or (y0>self.yr[1]):
                    return [],[],[],[]
            else:
                pass

        # disperse each polygon vertix
        xg,yg=self.disperse(xd,yd,wav,band=band)


        # apply clipping
        xg=np.clip(xg,0,self.naxis[0])
        yg=np.clip(yg,0,self.naxis[1])

        # clip against a pixel grid
        x,y,area,indices=self.polyclip(xg,yg)

        # replicate the wavelength indices
        n=len(x)
        if n>0:
            i0=indices[0:-1]
            i1=indices[1:]
            gg=np.where(i1 != i0)[0]
            lam=np.empty(n,dtype=np.uint16)
            for g,a,b in zip(gg,i0[gg],i1[gg]):
                lam[a:b]=g
        else:
            lam=[]


        return x,y,lam,area


    def sensitivity(self,wav,**kwargs):
        ''' method to evaluate the sensitivity function '''
        return self.sens(wav,**kwargs)


    
    
                
