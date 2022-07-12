#import h5py
import numpy as np
from skimage import measure
import pandas as pd

from .hdf5columns import HDF5Columns
from . import attributes
from ...config import Config

class HDF5Table(HDF5Columns):
    COLORS={'0':'white','+1':'#1f77b4','+2':'#ff7f0e','+3':'#2ca02c',
            '+4':'#d62728','+5':'#9467bd','+6':'#8c564b','+7':'#e377c2',
            '+8':'#7f7f7f','+9':'#bcbd22','+10':'#17becf','-1':'#aec7e8',
            '-2':'#ffbb78','-3':'#98df8a','-4':'#ff9896','-5':'#c5b0d5',
            '-6':'#c49c94','-7':'#f7b6d2','-8':'#c7c7c7','-9':'#dbdb8d',
            '-10':'#9edae5'}    
        
    def __init__(self,dims=None,**kwargs):
        self.dims=dims
        self.attrs={k:v for k,v in kwargs.items()}


        self.compargs=Config().h5pyargs
        
        for column in self.COLUMNS:
            self[column]=list()

        if self.dims:
            # this presumes that the dimensions are in the same order as
            # COLUMNS.  This would break if COLUMNS variable is reordered
            for name,dim in zip(self.COLUMNS,self.dims):
                self.attrs[f'n{name}']=self.DTYPES[name](dim)

    def __str__(self):
        return f'{self.__class__.__name__} for {self.name}'

    def __len__(self):
        return len(self[self.COLUMNS[0]])

    def __setitem__(self,k,v):
        super().__setitem__(k,list(v))


    def __mul__(self,a):
        self['val']=self.get('val')*a
        return self

    def __rmul__(self,a):
        return self.__mul__(a)

    def __imul__(self,a):
        return self.__mul__(a)

    def __itruediv__(self,a):
        self['val']=self.get('val')/a
        return self

    def __iter__(self):
        yield from zip(*self.values())



    @property
    def columns(self):
        return tuple(self.keys())

    @property
    def ncolumns(self):
        return len(self.keys())


    def clear(self):
        for v in self.values():
            v.clear()

    def extend(self,*args):
        for column,arg in zip(self.columns,args):
            self[column].extend(arg)

    def get(self,k,dtype=None):
        if dtype is None:
            dtype=self.DTYPES[k]
        return np.array(self[k],dtype=dtype)
                

    def wavelengths(self,lam=None):
        if 'wav0' in self.attrs and 'dwav' in self.attrs:
            wav0=self.attrs['wav0']
            dwav=self.attrs['dwav']

            if lam is None:
                lam=self.get('lam')

            wav=wav0+lam*dwav
            return wav



    def compute_xyg(self):
        LOGGER.debug('make this use utilities')

        xyg=np.ravel_multi_index((self.get('x'),self.get('y')),
                                 dims=(self.attrs['nx'],self.attrs['ny']),
                                 order='F')
        return xyg

    def as_pandas(self):
        data={column: self.get(column) for column in self.COLUMNS}
        data=pd.DataFrame(data=data)
        return data

    def as_array(self):
        # get the right dtypes
        dtype=[(column,self.DTYPES[column]) for column in self.COLUMNS]

        # create an array
        data=np.empty((len(self),),dtype=dtype)

        # fill the array
        for column in self.COLUMNS:
            data[column]=self.get(column)
        return data

    def bounding_box(self,dx=(0,0),dy=(0,0)):
        if 'x' in self and 'y' in self:
            x=self.get('x')
            x0=np.amin(x)-dx[0]
            x1=np.amax(x)+dx[1]
            if 'nx' in self.attrs:
                x0=max(x0,0)
                x1=min(x1,self.attrs['nx'])

            y=self.get('y')
            y0=np.amin(y)-dy[0]
            y1=np.amax(y)+dy[1]
            if 'ny' in self.attrs:
                y0=max(y0,0)
                y1=min(y1,self.attrs['ny'])
        else:
            x0=x1=y0=y1=None

        return x0,x1,y0,y1


    def select(self,g):
        for k in self.columns:
            self[k]=self.get(k)[g]

    def threshold(self,a):
        g=self.get('val')>=a
        self.select(g)


    def region(self,order=None,mask=False,size=12):
        if ('x' in self.COLUMNS) and ('y' in self.COLUMNS):
            x=self.get('x')
            y=self.get('y')

            # get range
            x0,x1=np.amin(x)-1,np.amax(x)+1
            y0,y1=np.amin(y)-1,np.amax(y)+1
            
            # make a binary image
            shape=(y1-y0+1,x1-x0+1)
            img=np.zeros(shape,dtype=float)        
            img[y-y0,x-x0]=1
            
            # contour the image
            contours = measure.find_contours(img,level=0.5)
            
            # reset the contours
            cy=contours[0][:,0]+y0
            cx=contours[0][:,1]+x0
            
            coord=','.join('{},{}'.format(*xy) for xy in zip(cx,cy))
            font=f'helvetica {int(size)}'


            color=self.COLORS.get(order,'black')
            

            reg=f'{int(mask)}polygon({coord}) # color={color} '+\
                f'text={{{self.name}}} edit=0 move=0 rotate=0 fixed=1 '+\
                f'font="{font}" width={int(width)}'

            #if family not in ('helvitica','times','courier'):
            #    family='helvetica'
            #font=f'{family} {int(size)}'

            #if bold:
            #    font+=' bold'
            #if italic:
            #    font+= 'italic'
            
            #reg=f'{int(mask)}polygon({coord}) # color={color} '+\
            #    f'text={{{self.name}}} edit={int(edit)} move={int(move)} '+\
            #    f'rotate={int(rotate)} width={int(width)} font="{font}" '+\
            #    f'fixed={int(fixed)}'


        else:
            LOGGER.error('Cannot make region, table does not "x" and "y"')
            
        return reg



    def write_hdf5(self,h5):        
        hd=h5.create_dataset(self.name,data=self.as_array(),**self.compargs)
        for k,v in self.attrs.items():
            attributes.write(hd,k,v)
