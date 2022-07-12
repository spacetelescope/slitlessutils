import numpy as np
from dataclasses import dataclass

from ...utilities import headers

@dataclass
class Grating:
    grating: str
    pupil: str
    units: str
    wave0: float
    wave1: float
    dwave: float


    @classmethod
    def from_xml(cls,xml,grating):
        kls=cls._get_class(grating)
        obj=kls(grating,xml.attrib['name'],xml.attrib['units'],
            float(xml.attrib['wave0']),float(xml.attrib['wave1']),float(xml.attrib['dwave']))
        return obj

    @classmethod
    def from_header(cls,hdr):
        filt=hdr['FILTER']

        kls=cls._get_class(filt)
        obj=kls(filt,hdr.get('PUPIL',''),hdr.get('WUNIT','A'),
            hdr['WAVE0'],hdr['WAVE1'],hdr['DWAVE'])

        return obj

    @staticmethod
    def _get_class(filt):
        f=filt[0].lower()
        if f=='g':
            return Grism
        elif f=='p':
            return Prism
        else:
            raise NotImplementedError(f'Grating format {filt} is unknown')

    def __bool__(self):
        return self.wave0 is not None and \
            self.wave1 is not None and \
            self.dwave is not None and \
            self.wave1 > self.wave0 and \
            self.wave0 > 0. and \
            self.wave1 > 0. and \
            self.dwave > 0.

    def update_header(self,hdr):
        hdr.set('GTYPE',value=self.GTYPE,comment='slitless grating type')
        hdr.set('wave0',value=self.wave0,comment=f'Initial wavelength')
        hdr.set('wave1',value=self.wave1,comment=f'Final wavelength')
        hdr.set('dwave',value=self.dwave,comment=f'Sampling wavelength')
        hdr.set('wunit',value=self.units,comment='units on wavelength')
        headers.add_stanza(hdr,'Wavelength Parameters',before='GTYPE')


class Grism(Grating):
    GTYPE = 'grism'

    def __call__(self,i,nsub=1):
        return self.wave0+i*self.dwave/nsub

    def __len__(self):
        return int(np.ceil((self.wave1-self.wave0)/self.dwave)) + 1

    def indices(self,wav):
        w=wav-self.wave0+self.dwave/2.
        ind=np.floor(w/self.dwave).astype(int)
        ind=np.clip(ind,0,len(self)-1)
        return ind


    def limits(self,nsub=1):
        w0=self.wave0-self.dwave/2.
        w1=self.wave1+self.dwave/2.
        n=int((w1-w0)*nsub/self.dwave)+1
        return np.linspace(w0,w1,n)

    def wavelengths(self,nsub=1):
        delta=self.dwave/nsub
        wav=np.arange(self.wave0,self.wave1+delta,delta,dtype=float)
        return wav



class Prism(Grating):
    GTYPE ='prism'

    ''' Specify a Prism spectral element.

    Key thing to know here is that the dispersion is geometric.  Could use the
    native numpy, but I want to specify the endpoints and the linear separation
    between the first two elements as dwave

    wave = wave0 *(wave1/wave0)**(b*i)

    where

    b=np.log(1+dwave/wave0)/np.log(wave1/wave0)

    and then the total number of possible elements will be

    n=1+int(np.floor(1/b))

    '''

    # a=self.wave1/self.wave0
    # wav = self.wave0 * a**(b*i)
    # b=np.log(1+self.dwave/self.wave0)/np.log(a)
    # n=1+int(np.floor(1/b))
    #
    # ind=np.floor(np.log(wav/self.wave0)/np.log(1+self.dwave/self.wave0)).astype(int)


    def __call__(self,i,nsub=1):
        pass



    # def __call__(self,i,nsub=1):
    #     n=int((self.wave1-self.wave0)/self.dwave)+1
    #     a=1+np.log(self.wave1/self.wave0)/n
    #     return self.wave0*a**i
    #
    # def __len__(self):
    #     w0=self.wave0-self.dwave/2.
    #     w1=self.wave1+self.dwave/2.
    #     return int((w1-w0)/self.dwave)+1

    def indices(self,wav):
        raise NotImplementedError


    # def limits(self,nsub):
    #     w0=self.wave0-self.dwave/2.
    #     w1=self.wave1+self.dwave/2.
    #     n=int((w1-w0)*nsub/self.dwave)+1
    #     return np.geomspace(w0,w1,n)

    def wavelengths(self,nsub=1):
        raise NotImplementedError


if __name__=='__main__':
    hdr={'FILTER':'G102','PUPIL':'none','WUNIT':'angstrom',
        'WAVE0':7500.,'WAVE1':12500.,'DWAVE':25.}

    g=Grating.from_header(hdr)
    print(g)
