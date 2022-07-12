from astropy.io import fits


from ....logger import LOGGER
from ...astrometry import WCS


class WFSSImage(WCS):
    TYPE = ''             # will be 'simulated' or 'observed'

    def __init__(self,hdr,detconf):
        self.hdr=hdr
        WCS.__init__(self,self.hdr)

        # copy some stuff over from the config object
        self.name=detconf.name
        self.bitmask=detconf.bitmask
        self.extensions={k:v for k,v in detconf.extensions.items()}



    @property
    def shape(self):
        return self.hdr['NAXIS2'],self.hdr['NAXIS1']

    @property
    def npix(self):
        return self.hdr['NAXIS1']*self.hdr['NAXIS2']

    def __str__(self):
        return f'WFSS image for {self.name}'


class WFSSFile(dict):

    GZIP = False
    TYPE = ''              # will be 'simulated' or 'observed'

    def __init__(self,dataset,insconf):

        # just copy some data over from the config object
        self.dataset=dataset
        self.telescope=insconf.telescope
        self.instrument=insconf.instrument
        self.grating=insconf.grating
        self.pupil=insconf.pupil
        self.bunit=insconf.bunit
        self.suffix=insconf.suffix

        
    def __str__(self):
        return f'{self.TYPE} WFSS image for {self.dataset}'

    @property
    def setting(self):
        return (self.telescope,self.instrument,self.grating,self.pupil)

    @property
    def filename(self):
        if self.GZIP:
            gzip='.gz'
        else:
            gzip=''

        if self.suffix=='':
            suffix=''
        else:
            suffix='_'+self.suffix
            
        return f'{self.dataset}{suffix}.fits{gzip}'
