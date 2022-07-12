import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d

from .hdf5table import HDF5Table

from ..utilities import indices

class OPT(HDF5Table):

    ''' class to hold an Object Profile Table (OPT) '''
    
    COLUMNS=('x','y','wav','val')


    def __init__(self,segid,dims=None,**kwargs):
        HDF5Table.__init__(self,dims=dims,**kwargs)
        
        self.segid=segid
        self.pdts={}
        self.pixels=[]


        

        
    @property
    def name(self):
        return str(self.segid)

    def append(self,pdt):
        pixel=pdt.pixel
        self.pdts[pixel]=pdt
        self.pixels.append(pixel)

    def decimate(self):
        if self.pdts:
            
            x,y=[],[]
            wav,val=[],[]
            
            for pdt in self.pdts.values():
                x.extend(pdt['x'])
                y.extend(pdt['y'])
                wav.extend(pdt.wavelengths())
                val.extend(pdt['val'])



            if len(x)>0:
                # ok... let's just save some space
                self.pdts.clear()

                # chagne datatypes
                x=np.array(x,dtype=int)
                y=np.array(y,dtype=int)
                wav=np.array(wav,dtype=int)
                val=np.array(val,dtype=float)
                
            
                # do the summations
                prof,x1,y1=indices.decimate(val,x,y,dims=self.dims)
                wave,x2,y2=indices.decimate(val*wav,x,y,dims=self.dims)
                wave/=prof

                
                if np.array_equal(x1,x2) and np.array_equal(y1,y2):
                    # stuff it back into the OPT
                    self.clear()
                    self['x']=x1
                    self['y']=y1
                    self['wav']=wave
                    self['val']=prof

                else:
                    raise RuntimeError("invalid x/y arrays")
            else:
                pass



    def as_image(self,show=True):
        x=self.get('x')
        y=self.get('y')
        
        x0,x1=np.amin(x),np.amax(x)
        y0,y1=np.amin(y),np.amax(y)

        dx=x1-x0+1
        dy=y1-y0+1
        
        
        prof=np.full((dy,dx),np.nan,dtype=float)
        wave=np.full((dy,dx),np.nan,dtype=float)
        
        x-=x0
        y-=y0

                
        prof[y,x]=self.get('val')
        wave[y,x]=self.get('wav')
        
   
        #xx=np.arange(dx)
        #yy=np.arange(dy)
        #W = interp2d(x, y, wave[y,x],fill_value='extrapolate')
        #wave=W(xx,yy)
        
        

        if show:
            padding=((3,3),(3,3))
            wave=np.pad(wave,padding,mode='constant',constant_values=np.nan)
            prof=np.pad(prof,padding,mode='constant',constant_values=np.nan)
        
        
            w0,w1=np.nanmin(wave),np.nanmax(wave)
            levels=np.linspace(w0,w1,10)
        
            fig,axes = plt.subplots(1,1,figsize=(10,3))
            axes.imshow(prof,origin='lower',interpolation='nearest')

            
            #,aspect=1)
            axes.contour(wave,levels,cmap='Spectral_r',origin='lower')

            #plt.tight_layout()
            fig.canvas.draw()


            #xlab=axes.get_xticklabels()
            #for i in range(len(xlab)):
            #    try:
            #        v=int(xlab[i].get_text())
            #    except:
            #        v=-1*int(xlab[i].get_text()[1:])
            #    xlab[i].set_text(str(v+x0))
            #axes.set_xticklabels(xlab)
            #            
            #ylab=axes.get_yticklabels()
            #for i in range(len(ylab)):
            #    try:
            #        v=int(ylab[i].get_text())
            #    except:
            #        v=-1*int(ylab[i].get_text()[1:])
            #    ylab[i].set_text(str(v+y0))
            #axes.set_yticklabels(ylab)
                        

            plt.show()
            

        
        return prof,wave
        
        
        
        
        
