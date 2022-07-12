from astropy.io import fits
import pandas as pd
import numpy as np

from ..photometry import Throughput
from ..utilities import headers
from ...logger import LOGGER

class DirectImages:
    ''' Class to contain the direct image properties '''

    def __init__(self,detindex=0):
        self.detindex=detindex
        self.filename=None

    def __len__(self):
        return len(self._data)

    def __iter__(self):
        for index,row in self._data.iterrows():
            yield tuple(row)

    def __str__(self):
        return f'Direct Images: \n {self._data}'

    @property
    def detband(self):
        return self._data.loc[self.detindex]['throughput']

    @property
    def detname(self):
        return self._data.loc[self.detindex]['filtfile']

    @property
    def detzero(self):
        #return self.detband.zeropoint
        return self.detband.zeropoint

    @classmethod
    def from_list(cls,filenames,**kwargs):
        obj=cls(**kwargs)

        # read the bandpasses and get the pivot wavelengths
        data={'filename':filenames,
              'photplam':[],
              'filtfile':[],
              'throughput':[]}
        for filename in data['filename']:
            hdr=fits.getheader(filename)
            if 'filtfile' in hdr:
                filtfile=hdr['filtfile']
            else:
                telescope=hdr['telescop']
                instrument=hdr['instrume']
                bandname=hdr['FILTER']

                filtfile=f'{telescope}_{instrument}_{bandname}.fits'

            band=Throughput.from_fits(filtfile)
            data['photplam'].append(band.photplam)
            data['throughput'].append(band)
            data['filtfile'].append(filtfile)

        # put the data in a pandas dataframe and sort
        obj._data=pd.DataFrame.from_dict(data)
        obj._data.sort_values('photplam',inplace=True)


        return obj


    @classmethod
    def from_file(cls,filename,**kwargs):
        files=np.loadtxt(filename,usecols=(0,),unpack=True)
        obj=cls().from_list(files,**kwargs)
        obj.filename=filename
        return obj

    def update_header(self,hdr):
        hdr['OBSFILE']=(self.filename,"Name of direc images' file")
        hdr['DETFILT']=(self.detname,'Name of detection filter')
        hdr['DETPLAM']=(self.detband.photplam,'pivot wavelength of detection band')
        hdr['DETZERO']=(self.detband.zeropoint,'AB mag zeropoint for detection band')
        hdr['NDIRECT']=(len(self),'Number of direct images')
        headers.add_stanza(hdr,'Direct Images',before='OBSFILE')
