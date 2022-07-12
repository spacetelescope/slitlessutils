import numpy as np
from astropy.io import fits
from contextlib import nullcontext


from ...logger import LOGGER
from ..photometry import Throughput
from ..astrometry import AstroImage
from ..utilities import indices,headers
from .source import Source
from .sedfile import SEDFile

class SourceCollection(dict):
    def __init__(self,segfile,detfile,maglim=np.inf,sedfile=None,**kwargs):
        LOGGER.info(f'Loading sources from {segfile} and {detfile}')

        self.segfile=segfile
        self.detfile=detfile
        self.sedfile=sedfile

        # get a magnitude limit
        self.maglim=maglim


        # what context manager to use
        if isinstance(self.sedfile,str):
            sedcontext=SEDFile(self.sedfile)
        else:
            sedcontext=nullcontext()


        # open the fits files and sort out their nature
        fitsargs={'mode':'readonly','ignore_missing_simple':True}
        with fits.open(self.segfile,**fitsargs) as hdus,\
            fits.open(self.detfile,**fitsargs) as hdud,\
            (SEDFile(self.sedfile) if isinstance(self.sedfile,str) else
                nullcontext()) as sedlib:

            # sort out the input nature
            nhdus,nhdud=len(hdus),len(hdud)
            if nhdus != nhdud:
                msg=f'Segmentation {self.segfile} and direct '+\
                    f'{self.obsdata.detfile} have differing number of '+\
                    f'extensions: {nhdus},{nhdud}'
                LOGGER.error(msg)
                return

            # is it an MEF?
            self.mef = nhdus > 1
            if self.mef:
                self._from_mef(hdus,hdud,sedlib,**kwargs)
            else:
                self._from_classic(hdus,hdud,sedlib,exten=0,**kwargs)


    def _from_classic(self,hdus,hdud,sedlib,exten=0,**kwargs):
        LOGGER.info(f'Loading classic segmentation map: {self.segfile}')

        # load image and segmap
        img=AstroImage.from_HDU(hdud[exten])
        seg=AstroImage.from_HDU(hdus[exten])


        # find pixels for each object
        ri=indices.reverse(seg.image,ignore=(0,))
        for segid,(y,x) in ri.items():

            # get a bounding Box
            x0,x1=np.amin(x),np.amax(x)
            y0,y1=np.amin(y),np.amax(y)

            # extract the subimages
            subimg=img.extract(x0,x1,y0,y1)
            subseg=seg.extract(x0,x1,y0,y1)

            # load the source
            source=Source(subimg,subseg,segid,zeropoint=img.zeropoint)

            if source and source.mag < self.maglim:
                self[segid]=source
                if sedlib:
                    source.load_seds(sedlib,throughput=img.throughput)


    def set_extraction_parameters(self):
        pass



    def update_header(self,hdr):
        maglim='INF' if np.isinf(self.maglim) else self.maglim

        sedfile=self.sedfile if self.sedfile else ' '*8

        hdr.set('SEGFILE',value=self.segfile,comment='segmentation image')
        hdr.set('DETFILE',value=self.detfile,comment='direct image for profile weights')
        hdr.set('MEF',value=self.mef,comment='is it a MEF file?')
        hdr.set('SEDFILE',value=sedfile,comment='file *INPUT* containing SEDs')
        hdr.set('NSOURCE',value=len(self),comment='number of sources')
        hdr.set('MAGLIM',value=maglim,comment='magnitude limit of sources')
        headers.add_stanza(hdr,'Source Properties',before='SEGFILE')



    # def load_seds(self,sedfile,**kwargs):
    #     LOGGER.debug('normalization is not done!')
    #     with SEDFile(self.sedfile) as sedlib:
    #         for source in self.values():
    #             source.load_seds(sedlib,**kwargs)

    def write_seds(self,**kwargs):
        for source in self.values():
            source.write_seds(**kwargs)


    def __str__(self):
        return f'Source Collection with {len(self)} sources'

        #indent=' '
        #
        #cont=lambda i,n: '\u2515' if i==n-1 else '\u251D'
        #left=lambda i,n: ' ' if i==n-1 else '\u2502'
        #
        #s=f'Source Collection:'
        #
        #n=len(self):
        #for i,(segid,source) in enumerate(self.items()):
        #
        #    s+=f'\n{indent}{cont(i,n)} {segid}'
