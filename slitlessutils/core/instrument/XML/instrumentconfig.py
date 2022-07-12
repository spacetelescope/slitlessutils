
import xml.etree.ElementTree as ET
from astropy.io import fits
import numpy as np
import os


from .detectorconfig import DetectorConfig
from .siaf import SIAF
from .grating import Grating
from ...utilities import headers
from ....config import Config


class InstrumentConfig(dict):

    def __init__(self,telescope,instrument,grating,xmlfile=None,
                 pupil='',detectors=None,orders=None):

        if xmlfile is None:
            xmlfile='instruments.xml'

        self.xmlfile=Config().get_reffile(xmlfile,path='instruments')
        tree=ET.parse(self.xmlfile)
        xml=tree.getroot()


        tel=self.get_node(xml,telescope)
        if tel is not None:
            ins=self.get_node(tel,instrument)
            if ins is not None:
                self.bunit=ins.attrib.get('bunit','ELECTRONS/S')
                self.suffix=ins.attrib.get('suffix','flt')

                grt=self.get_node(ins,grating)

                if grt is not None:
                    pup=self.get_node(grt,pupil)
                    if pup is not None:
                        self.telescope=telescope
                        self.instrument=instrument
                        self.grating=grating
                        self.pupil=pupil

                        self.refsiaf=SIAF.from_xml(ins)
                        self.parameters=Grating.from_xml(pup,grating)

                        self.header=fits.Header()
                        for hdr in ins.findall('header'):
                            self.header[hdr.attrib['keyword']]=(hdr.text,hdr.attrib['comment'])

                        # get all the detector names if not specified
                        for det in ins.findall('detector'):
                            name=det.attrib['name']

                            if detectors is None:
                                self[name]=DetectorConfig.from_xml(det,self.refsiaf,grating,pupil,orders=orders)
                            else:
                                if name in detectors:
                                    self[name]=DetectorConfig.from_xml(det,self.refsiaf,grating,pupil,orders=orders)


        self.orders = self[list(self.keys())[0]].orders



    def update_header(self,hdr):
        ''' method to update a fits header '''



        # record the telescope/instrument combination
        hdr['TELESCOP']=(self.telescope,'telescope used to acquire the data')


        # do something different for each insturment
        if self.instrument in ("WFC3IR",'WFC3UVIS'):
            #hdr['INSTRUME']=('WFC3','identifier for instrument used to acquire data')
            hdr['FILTER']=(self.grating,'element selected from filter wheel')
        elif self.instrument in ("ACSWFC","ACSSBC"):
            #hdr['INSTRUME']=("ACS",'identifier for instrument used to acquire data')
            hdr['FILTER1']=(self.grating,'element selected from filter wheel 1')
            hdr['FILTER2']=('CLEAR2L','element selected from filter wheel 2')
        else:
            hdr['INSTRUME']=(self.instrument,'identifier for instrument used to acquire data')
            hdr['FILTER']=(self.grating,'element selected from filter wheel')


        # put in some things for observing mode
        hdr['OBSTYPE']=('SPECTROSCOPIC','observation type - imaging or spectroscopic')
        hdr['PUPIL']=(self.pupil.ljust(8,' '),'pupil filter')

        # let's add some extra ones that came with the config.
        for k in self.header.keys():
            hdr[k]=(self.header[k],self.header.comments[k])


        headers.add_stanza(hdr,'INSTRUMENT CONFIGURATION INFORMATION',before='TELESCOP')



    @property
    def name(self):
        return self.instrument


    @staticmethod
    def get_node(xml,name):
        for x in xml:
            if x.attrib.get('name')==name:
                return x



if __name__=='__main__':
    ins=InstrumentConfig('HST','WFC3UVIS','G280')
