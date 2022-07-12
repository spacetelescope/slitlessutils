from astropy.io import fits
from datetime import datetime
import numpy as np

from ...utilities import headers
from .wfss import WFSSFile,WFSSImage
from ....logger import LOGGER

class SimulatedImage(WFSSImage):
    ''' Class to hold a Simulated WFSS Image '''

    TYPE='simulated'

    def __init__(self,hdr,detconf,exptime,background=0.0,bunit='ELECTRONS/S'):
        WFSSImage.__init__(self,hdr,detconf)


        self.noise=detconf.noise
        self.exptime=exptime
        self.background=background
        self.bunit=bunit



    def make_header(self,imgtype):
        ''' make a default header '''
        exten=self.extensions[imgtype]
        hdr=super().mkhdr(exten.dtype)

        exten.update_header(hdr)

        return hdr

    def make_HDUs(self,sci,addnoise=True):#,noisepars):
        ''' make the HDUs for a simulated image '''


        if addnoise:
            # add in Poisson noise terms
            pvars=sci*self.exptime
            pvarc=(self.noise.dark+self.background)*self.exptime
            sci=np.random.poisson(lam=pvars+pvarc,size=sci.shape)-pvarc

            # add in Gaussian noise terms
            gvar=self.noise.read*self.noise.read
            sci=np.random.normal(loc=sci,scale=np.sqrt(gvar))/self.exptime

            # create an uncertainty image
            unc=np.sqrt(pvars+pvarc+gvar)/self.exptime
        else:
            unc=np.ones_like(sci)


        # change the units to match BUNIT
        bunit=self.bunit.lower()
        if bunit in ('electron','electrons','e','e-'):
            sci *= self.exptime
            unc *= self.exptime

        # make some HDUs
        scihdu=fits.ImageHDU(data=sci.astype(self.extensions['science'].dtype),
                             header=self.make_header('science'))
        unchdu=fits.ImageHDU(data=unc.astype(self.extensions['uncertainty'].dtype),
                             header=self.make_header('uncertainty'))
        dqahdu=fits.ImageHDU(data=np.ones_like(sci,dtype=self.extensions['dataquality'].dtype),
                             header=self.make_header('dataquality'))

        # add the BUNIT --- this is hardcoded for HST
        scihdu.header['BUNIT']=(self.bunit,'brightness unit')
        unchdu.header['BUNIT']=(self.bunit,'brightness unit')
        dqahdu.header['BUNIT']=('UNITLESS','brightness unit')




        # update some headers
        self.noise.update_header(scihdu.header,addnoise=addnoise)


        return scihdu,unchdu,dqahdu

class SimulatedFile(WFSSFile):
    ''' A Class to hold a Simulated WFSS File '''

    TYPE='simulated'

    def __init__(self,dataset,crval1,crval2,orientat,exptime,
                 insconf,targname='',**kwargs):
        WFSSFile.__init__(self,dataset,insconf)

        # make a primary hdr
        self.phdr=fits.Header()

        # update the primary header
        insconf.update_header(self.phdr)

        # put in some things about the targname
        self.phdr['TARGNAME']=(targname.ljust(8,' '),"proposer's target name")
        headers.add_stanza(self.phdr,'Target Information',before='TARGNAME')


        # add some exposure Information
        now=datetime.now()
        t0=2444493.50000-2_400_000.5
        self.phdr['DATE-OBS']=(now.strftime("%Y-%m-%d"),
                               'UT date of start of observation (yyyy-mm-dd)')
        self.phdr['TIME-OBS']=(now.strftime("%H:%M:%S"),
                               'UT time of start of observation (hh:mm:ss)')
        self.phdr['EXPSTART']=(t0,'exposure start time (Modified Julian Date)')
        self.phdr['EXPEND']=(t0+exptime/86400,'exposure end time (Modified Julian Date)')
        self.phdr['EXPTIME']=(exptime,'exposure duration (seconds)')
        self.phdr['EXPFLAG']=('NORMAL','Exposure interruption indicator')
        headers.add_stanza(self.phdr,'EXPOSURE INFORMATION',before='DATE-OBS')


        # now process each detector
        for detname,detconf in insconf.items():
            hdr=detconf.mkwcs(crval1,crval2,orientat)
            self[detname]=SimulatedImage(hdr,detconf,exptime,bunit=self.bunit,**kwargs)
