import numpy as np
from dataclasses import dataclass
from ...utilities import headers

@dataclass
class Noise:
    dark: float
    read: float
    #distribution: str

    @classmethod
    def from_xml(cls,xml):
        return cls(float(xml.attrib.get('dark',0.)),
                   float(xml.attrib.get('read',0.)))
    #xml.attrib.get('distribution','poisson'))

    def update_header(self,hdr,addnoise=True):
        hdr['NOISE']=(addnoise,'Was noise added?')
        if addnoise:
            hdr['RDNOISE']=(self.read,'assumed readnoise in e-')
            hdr['DARKRATE']=(self.dark,'assumed dark rate in e-/s')
            #hdr['UNCDIST']=(self.distribution,'uncertainty distribution model')
        headers.add_stanza(hdr,'Noise Properties',before='NOISE')

    # def add_noise(self,scihdu,unchdu,exptime,background,addnoise=True):
    #     LOGGER.debug('this is to be deprecated')
    #     if addnoise:
    #         # save the dtypes
    #         scitype,unctype=scihdu.dtype,unchdu.dtype
    #
    #
    #         # add in poisson terms
    #         pvars=scihdu.data*exptime
    #         pvarc=(background+self.dark)*exptime
    #         scihdu.data=np.random.poisson(pvars+pvarc)-pvarc
    #
    #         # add in the Gaussian terms
    #         gvar=self.read*self.read
    #         scihdu.data+= np.random.normal(scale=np.sqrt(gvar))
    #         scihdu.data=scihdu.data.astype(scitype)
    #
    #         # now update the uncertainty
    #         unchdu.data=(np.sqrt(pvars+pvarc+gvar)/exptime).astype(unctype)
    #
    #
    #             #if self.distribution=='poisson':
    #             #    scihdu.data=np.random.poisson(var)-self.dark*exptime
    #             ##elif self.distribution in ('normal','gaussian'):
    #             #scihdu.data=np.random.normal(scale=np.sqrt(var))
    #
    #
    #
    #             #var=(scihdu.data+self.dark)*exptime+(self.read*self.read)
    # 
    #             #unchdu.data=np.sqrt(var)/exptime
    #             #if isinstance(dtype,type):
    #             #    unchdu.data=unchdu.data.astype(dtype)
    #
    #             #if self.distribution=='poisson':
    #             #    scihdu.data+=(np.random.poisson(var)-var)/exptime
    #             #elif self.distribution in ('normal','gaussian'):
    #             #    scihdu.data+=np.random.normal(scale=unc.data)
    #
    #         # reassign the dtypes
    #         #scihdu.data=scihdu.data.astype(scitype)
    #         #unchdu.data=unchdu.data.astype(unctype)
    #     #else:
    #     #        self.update_header(scihdu.header,noise=False)
    #     #        self.update_header(unchdu.header,noise=False)
    #     else:
    #         self.update_header(scihdu.header,noise=False)
    #
