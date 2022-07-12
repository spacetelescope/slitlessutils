from astropy.io import fits
import pandas as pd
import os
from dataclasses import dataclass
from io import StringIO

from ...instrument import InstrumentConfig
from .observed import ObservedFile
from .simulated import SimulatedFile
from ....logger import LOGGER
from ...utilities import headers


@dataclass(frozen=True,eq=True)
class ConfKey:
    ''' a dataclass to hold configuration data '''
    telescope: str
    instrument: str
    grating: str
    pupil: str

    @classmethod
    def from_fits(cls,filename):
        # read primary header
        h0=fits.getheader(filename,ext=0)

        # name of the instrument
        inst=h0['INSTRUME']

        # Gotta worry about each instrument specially
        if inst == 'WFC3':
            inst+=h0['DETECTOR']
            filt=h0['FILTER']
        elif inst=='ACS':
            inst+=h0['DETECTOR']
            filt=h0['FILTER1']
        else:
            filt=h0['FILTER']


        # cross filter
        pupil=h0.get('PUPIL',None)

        # inialize the object
        return cls(h0['TELESCOP'],inst,filt,pupil)


    def __str__(self):
        s=f'{self.telescope}/{self.instrument} {self.grating}'
        if self.pupil is not None and self.pupil not in ('',' '):
            s+=f' +{self.pupil}'
        return s

    def load_config(self,**kwargs):
        insconf=InstrumentConfig(self.telescope,self.instrument,self.grating,\
            pupil=self.pupil,**kwargs)
        return insconf


class DataKey:
    ''' A base class to configure an observation '''
    pass

@dataclass
class ObservedKey(DataKey):
    ''' A dataclass to hold configurations for an ObservedFile '''
    filename: str

    def load_file(self,insconf):
        return ObservedFile(self.filename,insconf)

@dataclass
class SimulatedKey(DataKey):
    ''' A dataclass to hold configurations for a SimulatedFile '''
    dataset: str
    ra: float
    dec: float
    orientat: float
    exptime: float
    background: float

    def load_file(self,inscnf):
        return SimulatedFile(self.dataset,self.ra,self.dec,self.orientat,self.exptime,
                             inscnf,background=self.background)

################################################################################
################################################################################
################################################################################



class WFSSCollection(dict):
    ''' A collection of WFSS data (simualted or observed) '''

    COMMENTS = ('%','#','&','$','!',';')      # comments in the data

    def __init__(self,filename=None,filetype=None,orders=None):
            self.filename=filename
            self.filetype=filetype
            self.orders=orders

    @property
    def nconfig(self):
        return len(self)

    @property
    def nfiles(self):
        return sum(len(v) for v in self.values())

    @property
    def npixels(self):
        ''' return the total number of pixels in this collection '''
        npix=0
        for cnfkey,files in self.items():
            insconf=cnfkey.load_config(orders=[])  #self.orders)
            npix+=(insconf.npixels*len(files))
        return npix

    @classmethod
    def from_list(cls,filenames,**kwargs):
        ''' classmethod to load WFSS files from a python list of filenames '''

        LOGGER.info(f'Loading observed WFSS files from python list')

        filename=kwargs.get('filename','<LIST>')
        obj=cls(filename=filename,filetype='observed',**kwargs)
        for filename in filenames:
            if filename and os.path.exists(filename):
                cnfkey=ConfKey.from_fits(filename)
                datkey=ObservedKey(filename)
                obj[cnfkey]=datkey

        return obj


    @classmethod
    def from_file(cls,obsfile,**kwargs):
        ''' class method to load WFSS files from a file list '''
        LOGGER.info(f'Loading observed WFSS files from {obsfile}')

        with open(obsfile,'r') as fp:
            filesnames=[line.strip() for line in fp if line[0] not in cls.COMMENTS]

        obj=cls.from_list(filenames,filename=obsfile,**kwargs)

        return obj


    @classmethod
    def from_dataframe(cls,df,**kwargs):
        ''' classmethod to load WFSS files from a pandas DataFrame '''

        LOGGER.info(f'Loading simulated WFSS files from pandas.DataFrame')

        obj=cls(filetype='simulated',**kwargs)
        for r in df.itertuples():
            cnfkey=ConfKey(r.telescope,r.instrument,r.grating,r.pupil)
            datkey=SimulatedKey(r.Index,r.ra,r.dec,r.orientat,r.exptime,r.background)
            obj[cnfkey]=datkey
        return obj

    @classmethod
    def from_wcsfile(cls,wcsfile,**kwargs):
        ''' classmethod to load WFSS files from a CSV that contains WCS info '''

        LOGGER.info(f'Loading simulated WFSS images from {wcsfile}')

        # set some defaults (or constants)
        defaults={'telescope': kwargs.get('telescope','HST'),
                  'instrument': kwargs.get('instrument','WFC3IR'),
                  'grating': kwargs.get('grating','G102'),
                  'pupil': kwargs.get('pupil',''),
                  'ra': kwargs.get('ra',42.),                       # in deg
                  'dec': kwargs.get('dec',42.),                     # in deg
                  'orientat': kwargs.get('orientat',0.0),           # in deg
                  'exptime': kwargs.get('exptime',1000.),           # in sec
                  'background': kwargs.get('background',0.0)}       # in e/s

        # parse the CSV file, but do it this way to parse
        # constant settings set behind a comment card
        with open(wcsfile,'r') as fp:
            data=''
            for line in fp:
                if line[0] in cls.COMMENTS:
                    line=line[1:].strip()
                    tokens=[t.strip() for t in line.split('=')]
                    if len(tokens)==2:
                        defaults[tokens[0]]=tokens[1]
                else:
                    data+=line

        # now read it as a dataframe
        converters={k:type(v) for k,v in defaults.items()}
        df=pd.read_csv(StringIO(data),comment='#',converters=converters,
                       index_col='dataset')

        # add the constants back
        for k,v in defaults.items():
            if k not in df:
                df.insert(0,k,v,allow_duplicates=False)

        # now parse it as a dataframe
        return cls.from_dataframe(df,filename=wcsfile,**kwargs)


    def files(self):
        ''' an iterator to loop over files '''
        for cnfkey,datkeys in self.items():
            insconf=cnfkey.load_config(orders=self.orders)
            for datkey in datkeys:
                yield datkey.load_file(insconf)

    def __iter__(self):
        ''' an iterator to get useful info '''
        for cnfkey,datkeys in self.items():
            insconf=cnfkey.load_config(orders=self.orders)
            for datkey in datkeys:
                yield (insconf,datkey.load_file(insconf))

    def __next__(self):
        ''' implement next() '''
        if not hasattr(self,'_iter'):
            self._iter=iter(self)
        return next(self._iter)

    def __str__(self):
        ''' override the str printing format '''

        indent=' '

        cont=lambda i,n: '\u2515' if i==n-1 else '\u251D'
        left=lambda i,n: ' ' if i==n-1 else '\u2502'

        s=f'{self.filetype.capitalize()} WFSS Collection:'

        n=len(self)
        for i,(cnfkey,data) in enumerate(self.items()):

            s+=f'\n{indent}{cont(i,n)} {cnfkey}'
            m=len(data)
            for j,name in enumerate(data):
                s+=f'\n{indent}{left(i,n)} {cont(j,m)} {name}'
        return s

    def __setitem__(self,cnfkey,datkey):
        ''' overriding the default nature of dict() '''

        if isinstance(cnfkey,ConfKey) and isinstance(datkey,DataKey):
            if cnfkey not in self:
                super().__setitem__(cnfkey,[])
            self[cnfkey].append(datkey)

    def get_parameters(self):
        ''' return the default extraction parameters for this collection '''

        LOGGER.debug('not the best way to set a default extraction')

        cnfkeys=list(self.keys())
        insconf=cnfkeys[0].load_config(orders=[])
        return insconf.parameters


    def update_header(self,hdr):
        ''' update fits header with some useful information '''
        filename=self.filename if self.filename else ' '*8


        cnfkeys=list(self.keys())
        if len(cnfkeys)==1:
            inscnf=cnfkeys[0].load_config(orders=[])
            inscnf.update_header(hdr)
        else:
            raise NotImplementedError("update for multiple gratings")

        hdr.set('NWFSS',value=self.nfiles,comment='number of WFSS exposures')
        #hdr.set('INSFILE',value=cnf.XMLFILE,comment='instrument library')
        hdr.set('WFSSFILE',value=self.filename,comment='WFSS filename')
        headers.add_stanza(hdr,'WFSS Observations',before='NWFSS')
