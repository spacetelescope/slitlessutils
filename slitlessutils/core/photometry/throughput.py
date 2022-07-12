from astropy.io import fits
import os

from .band import Band
from ...logger import LOGGER
from ...config import Config

class Throughput(Band):
    def __init__(self,*args,**kwargs):
        Band.__init__(self,*args,**kwargs)


        
    @classmethod
    def from_keys(cls,telescope,instrument,band):
        LOGGER.info('loading throught from keys')
        filename=f'{telescope}_{instrument}_{band}.fits'.lower()
        filename=Config().get_reffile(filename,path='bandpasses')

        data,header=fits.getdata(filename,exten=1,header=True)

        obj=cls(data['wavelength'],data['transmission'],
                 unit=header.get('TUNIT1',''))

        obj.telescope=telescope
        obj.instrument=instrument
        obj.band=band
        obj.filename=filename

        return obj

    @classmethod
    def from_file(cls,filename):
        tokens=os.path.splitext(filename)
        ext=tokens[-1][1:]
        if ext == 'fits':
            data,header=fits.getdata(filename,exten=1,header=True)
            obj=cls(data['wavelength'],data['transmission'],
                    unit=header.get('TUNIT1',''))
            
        elif ext in ('dat','txt','filt','ascii'):
            LOGGER.debug("read ascii filt file")
            
        else:
            LOGGER.warning(f"File extension {ext} is unknown")
            return

        obj.filename=filename

        return obj



    def __str__(self):
        return f'Throughput curve for {self.telescope}/{self.instrument} {self.band}'
