from numpy as np
from astropy.io import fits


from ...logger import LOGGER
from ..astrometry import AstroImage
from ..utilities import indices,headers



class SourceCollection(dict):
    def __init__(self,segfile,detfile,maglim=None,sedfile=None,**kwargs):
        LOGGER.info(f'Loading sources from {segfile}')

        self.segfile=segfile

        # get a magnitude limit
        if maglim is None:
            maglim=np.inf
        self.maglim=maglim

        LOGGER.debug("SORT OUT DETFILE")
        self.detfile=detfile


        
        # open the fits files and sort out their nature
        mode='readonly'
        with fits.open(self.segfile,mode=mode) as hdus,\
             fits.open(self.detfile,mode=mode) as hdui:

            nhdus,nhdui=len(hdus),len(hdui)
            if nhdus != nhdui:
                msg=f'Segmentation {self.segfile} and direct '+\
                    f'{self.obsdata.detfile} have differing number of '+\
                    f'extensions: {nhdus},{nhdui}'
                LOGGER.error(msg)
                return

            # is it an MEF?
            self.mef = nhdus > 1
            if self.mef:
                self._from_mef(hdus,hdui,**kwargs)
            else:
                self._from_classic(hdus,hdui,**kwargs)

        # check some things
        if len(self)>0:

            # load some SEDs?
            self.sedfile=sedfile
            if self.sedfile is not None:
                self.load_seds(self.sedfile)


    def _from_classic(self,hdus,hdui,exten=0,**kwargs):
        LOGGER.info(f'Loading classic segmentation map: {self.segfile}')

        # load image and segmap
        img=AstroImage.from_HDU(hdui[exten])
        seg=AstroImage.from_HDU(hdus[exten])

        # get the zeropoint
        zeropoint=self.get_zeropoint(img)

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
            source=Source(subimg,subseg,segid,zeropoint=zeropoint)

            if source and source.mag < self.maglim:
                self[segid]=source
                


    def get_zeropoint(self,img):
        ''' extract the zeropoint from the file '''

        # if BUNIT is not specified, then assume e-/s
        bunit=img.get('BUNIT','ELECTRONS/S').lower()
        if bunit=='mjy/sr':
            # convert from MJy/sr to Jy/pix and divid by reference flux
            z=-2.5*np.log10((1e6/3631.)*(np.pi/180)**2*(img.pixelscale/3600)**2)
        elif bunit in ('electrons/s','electron/s','e/s','e-/s'):
            z=self.images.detzero
        else:
            z=0.0

        return z


    def update_header(self,hdr):
        hdr.set('NSOURCE',value=len(self),comment='number of sources')
        hdr.set('SEGFILE',value=self.segfile,comment='segmentation image')
        hdr.set('MAGLIM',value=self.maglim,comment='magnitude limit of sources')
        hdr.set('MEF',value=self.mef,comment='is it a MEF file?')
        hdr.set('SEDFILE',value=self.sedfile,comment='file containing SEDs')
        #self.obsdata.update_header(hdr)
        headers.add_stanza(hdr,'Source Properties',before='NSOURCE')

            

    def load_seds(self,sedfile):
        with SEDFile(self.sedfile) as sedlib:
            for source in self.values():
                source.load_seds(sedlib)

    def write_seds(self,**kwargs):
        for source in self.values():
            source.write_seds(**kwargs)

