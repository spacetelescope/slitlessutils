import numpy as np
from astropy.io import fits
import os

from ....logger import LOGGER
from ....config import Config
from ..wfssconfig import AsciiConfig,Order,load_flatfield
from .extension import Extension
from .noise import Noise
from .siaf import SIAF
from .sip import SIP
from .polyclip import polyclip

class DetectorConfig(dict):


    @classmethod
    def from_xml(cls,xml,refsiaf,grating,pupil,orders=None):
        obj=cls()

        # get teh name of the detector
        obj.name=xml.attrib['name']

        # get SIAF info from instrument setup
        obj.siaf=SIAF.from_xml(xml)
        obj.refsiaf=refsiaf

        # get the module to do simulated noise
        obj.noise=Noise.from_xml(xml.find('noise'))

        # get the bitmask setting
        bitmask=xml.attrib.get('bitmask')
        if bitmask is not None:
            bitmask=int(bitmask)
        obj.bitmask=bitmask

        # grab any header keywords
        obj.header=fits.Header()
        for hdr in xml.findall('header'):
            obj.header[hdr.attrib['keyword']]=(hdr.text,hdr.attrib['comment'])



        # get some properites of the pixel grid/tangent point
        as_list=lambda *keys: [xml.attrib[n] for n in keys]
        obj.naxis=np.array(as_list('naxis1','naxis2'),dtype=int)
        obj.scale=np.array(as_list('scale1','scale2'),dtype=float)
        obj.crpix=np.array(as_list('crpix1','crpix2'),dtype=float)


        # get the file structure information
        obj.extensions={ext.text: Extension.from_xml(ext) for ext in
                        xml.findall('extension')}

        # load SIP information from XML
        obj.sip=SIP.from_xml(xml.find('sip'))


        # now read the WFSS config
        cfg=xml.find('config')
        conffile=None
        for fil in cfg.findall('file'):
            if fil.attrib['grating']==grating and fil.attrib['pupil']==pupil:
                conffile=fil.text
        if conffile is None:
            print("GRATING/PUPIL not found")
        path=os.path.join('instruments',cfg.attrib['path'])
        conffile=Config().get_reffile(conffile,path=path)
        conf=AsciiConfig(conffile)


        # load the order data from the wfssconfig object
        if orders is None:          # default of which orders to load
            orders=conf.orders
        for order in orders:
            obj[order]=Order.from_dict(conf,order)


        # get the name of the flat field.  NB: don't load it because
        # that can be a heavy data prodcut, we just want the name.
        # the loading of the flatfield will be handled by a method below
        obj.ffname=None
        if 'FFNAME' in conf:
            ffname=os.path.join(conf.confpath,conf['FFNAME'])
            if os.path.exists(ffname):
                obj.ffname=ffname

        return obj


    def drizzle(self,xd,yd,wav,order,band=None):
        xg,yg=self[order].disperse(xd,yd,wav,band=band)

        # clip the pixels to be in the range
        xg=np.clip(xg,0,self.naxis[0])
        yg=np.clip(yg,0,self.naxis[1])

        # clip against a pixel grid
        x,y,area,indices=polyclip.multi(xg,yg,self.naxis[0],self.naxis[1])

        # decompose the indices as wavelength indices
        n=len(x)
        if n>0:
            i0=indices[0:-1]
            i1=indices[1:]
            gg=np.where(i1 != i0)[0]
            lam=np.empty(n,dtype=np.uint16)
            for g,a,b in zip(gg,i0[gg],i1[gg]):
                lam[a:b]=g
        else:
            lam=[]

        return x,y,lam,area




    def load_flatfield(self,unity=False):
        ''' return the flat field '''

        if unity or (self.ffname is None):
            ff=load_flatfield()
        else:
            ff=load_flatfield(self.ffname)
        return ff

    def mkwcs(self,ra,dec,orientat):
        ''' return a header with WCS information give a crude pointing '''

        # compute the CRVALs based on the SIAF and pointing
        crvals=self.siaf.crvals(ra,dec,orientat,self.refsiaf)

        # compute a rotation matrix
        ang=np.radians(orientat)
        cs,sn=np.cos(ang),np.sin(ang)
        R=np.array([[cs,-sn],[sn,cs]])

        # compute a pixel matrix
        P=np.diag([-self.scale[0]/3600.,self.scale[1]/3600.])

        # compute CD matrix as dot product of the pixel and rotation matrices
        CD=np.dot(R,P)

        # make a fits header
        h=fits.Header()
        h['NAXIS']=(2,'number of axes')
        h['NAXIS1']=(self.naxis[0],'number of pixels in x')
        h['NAXIS2']=(self.naxis[1],'number of pixels in y')
        h['CRPIX1']=(self.crpix[0],'x-coordinate of reference pixel')
        h['CRPIX2']=(self.crpix[1],'y-coordinate of reference pixel')
        h['CRVAL1']=(crvals[0],'first axis value at reference pixel')
        h['CRVAL2']=(crvals[1],'second axis value at reference pixel')
        h['CDELT1']=(1.0,' ')
        h['CDELT2']=(1.0,' ')
        h['CD1_1']=(CD[0,0],'partial of first axis coordinate w.r.t. x')
        h['CD1_2']=(CD[1,0],'partial of first axis coordinate w.r.t. y')
        h['CD2_1']=(CD[0,1],'partial of second axis coordinate w.r.t. x')
        h['CD2_2']=(CD[1,1],'partial of second axis coordinate w.r.t. y')
        h['CTYPE1']=('RA---TAN','the coordinate type for the first axis')
        h['CTYPE2']=('DEC--TAN','the coordinate type for the second axis')
        h['EQUINOX']=(2000.,'equinox of coordinates')
        h['LATPOLE']=(90.,' ')
        h['LONGPOLE']=(180.,' ')
        h['ORIENTAT']=(-orientat,'position angle of image y axis (deg. e of n)')
        self.sip.update_header(h)

        return h


    @property
    def orders(self):
        ''' return the orders present as a list '''
        return list(self.keys())

    @property
    def npixels(self):
        ''' return the total number of pixels in this detector '''
        return self.naxis[0]*self.naxis[1]
