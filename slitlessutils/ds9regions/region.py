#from .font import Font
from .attributes import Attributes

class Region:
    DELIM=' #'
    def __init__(self,background=False,**kwargs):

        # initialize a region's attributes as all None
        self.attrs=Attributes.none()
        for k,v in kwargs.items():
            self.attrs[k]=v

        # region specific attributes
        self.attrs['background']=background
        self._unit=''

        
    def ds9format(self,attrs={},fk5=None):
        if isinstance(fk5,bool):
            self.fk5=fk5
        
        # format the attributes
        a=[]
        for k,v in self.attrs.items():
            if v is not None and v!=attrs.get(k,None):
                a.append(self.attrs.format_value(k))
        attrib=' '.join(filter(None,a))
        
        return f'{self.include}{self.region}{self.DELIM} {attrib}'
        
        
    def __str__(self):
        return self.ds9format()

    @property
    def include(self):
        if self.attrs['include'] is None or self.attrs['include']:
            return ''
        else:
            return'-'

        
    @property
    def unit(self):
        return self._unit
    
    
        
    @property
    def fk5(self):
        return self._fk5

    @fk5.setter
    def fk5(self,_fk5):
        if isinstance(_fk5,bool):
            self._fk5=_fk5
            self._unit='"' if self._fk5  else ''
        else:
            raise ValueError("fk5 must be a boolean")

