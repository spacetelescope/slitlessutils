
import numpy as np

import pysiaf
from dataclasses import dataclass
from ...utilities import headers


@dataclass
class SIAF:
    v2: float
    v3: float
    v3y: float
    naxis1: int
    naxis2: int

    @classmethod
    def from_xml(cls,xml):
        
        obj=cls(float(xml.attrib['v2']),
                float(xml.attrib['v3']),
                float(xml.attrib['v3y']),
                int(xml.attrib.get('naxis1',0)),
                int(xml.attrib.get('naxis2',0)))
        return obj

    def update_header(self,hdr,reference=False):
        ref='Reference' if reference else ''
        hdr.set('V2',value=self.v2,comment='V2 (arcsec)')
        hdr.set('V3',value=self.v3,comment='V3 (arcsec)')
        hdr.set('V3Y',value=self.v3y,comment='V3Y (deg)')
        headers.add_stanza(hdr,f'{ref} SIAF Info',before='V2')

    def __eq__(self,refsiaf):
        if isinstance(refsiaf,type(self)):
            return np.isclose(self.v2,refsiaf.v2) and \
                np.isclose(self.v3,refsiaf.v3) and \
                np.isclose(self.v3y,refsiaf.v3y)
        else:
            return False

    def crvals(self,ra,dec,orientat,refsiaf):

        # do a quick test, if the refsiaf == self, then
        # do not need to do anything
        if refsiaf == self:
            # if here, then the reference data is where we're at.  so,
            # let's just short-circute and use the ra/dec as given
            crvals=(ra,dec)
        else:
            # if here, then have to compute the CRVALS, so do this
            # based on Bryan Hilbert's email from Apr 15, 2020, I should
            # transform orientat to a local_roll, but that is really
            # complicated.  But the first order is:
            # local_roll=-orientat-self.v3y
            # is essentially that.  Higher order will require more coding.
            local_roll=-orientat-self.v3y
            A = pysiaf.rotations.attitude(refsiaf.v2, refsiaf.v3,
                                          ra,dec,local_roll)


            # compute the new positions
            crvals = pysiaf.utils.rotations.pointing(A, self.v2, self.v3)

        return crvals

        

        
@dataclass
class Grating:
    grating: str
    pupil: str
    wave0: float
    wave1: float
    dwave: float
    units: str

    @classmethod
    def from_xml(cls,xml,grating):
        obj=cls(grating,
                xml.attrib['name'],
                float(xml.attrib['wave0']),
                float(xml.attrib['wave1']),
                float(xml.attrib['dwave']),
                xml.attrib['units'])
        return obj

    @classmethod
    def from_header(cls,hdr,**kwargs):
        obj=cls(hdr['FILTER'],
                hdr.get('pupil',''),
                float(hdr['wave0']),
                float(hdr['wave1']),
                float(hdr['dwave']),
                hdr['wunit'])
        return obj

    #@_validate
    def __len__(self):
        return int(np.ceil((self.wave1-self.wave0)/self.dwave)) + 1


    #@_validate
    def limits(self,nsub=1):
        lim=np.arange(self.wave0-self.dwave/2,
                      self.wave1+self.dwave*(2+nsub)/(2*nsub),
                      self.dwave/nsub)
        return lim

    #@_validate
    def wavelengths(self,nsub=1):
        wav=np.arange(self.wave0,self.wave1+self.dwave/nsub,
                      self.dwave/nsub,dtype=float)
        return wav


    #@_validate
    def indices(self,nsub=1):
        ind=np.arange(len(self)*nsub+1,dtype=int)
        return ind


    def update_header(self,hdr):
        hdr.set('wave0',value=self.wave0,comment=f'Initial wavelength')
        hdr.set('wave1',value=self.wave1,comment=f'Final wavelength')
        hdr.set('dwave',value=self.dwave,comment=f'Sampling wavelength')
        hdr.set('wunit',value=self.units,comment='units on wavelength')
        headers.add_stanza(hdr,'Wavelength Parameters',before='wave0')


    def __call__(self,i,nsub=1):
        return self.wave0+i*self.dwav/nsub

    def __bool__(self):
        logic=self.wave0 is not None and \
            self.wave1 is not None and \
            self.dwave is not None and \
            self.wave1 > self.wave0 and \
            self.wave0 > 0. and \
            self.wave1 > 0. and \
            self.dwave > 0.
        return logic


