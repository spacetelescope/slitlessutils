from dataclasses import dataclass
import numpy as np
from ...utilities import headers

@dataclass
class Extension:
    name: str
    ver: int 
    dtype: type

    @classmethod
    def from_xml(cls,xml):
        try:
            dtype=eval(xml.attrib['dtype'])
        except:
            dtype=float
            
        return cls(xml.attrib['name'].upper(),
                   int(xml.attrib['ver']),
                   dtype)
                   

    
    def update_header(self,hdr):
        hdr['EXTNAME']=(self.name,'')
        hdr['EXTVER']=(self.ver,'')
        headers.add_stanza(hdr,f'Extension Data',before='EXTNAME')

