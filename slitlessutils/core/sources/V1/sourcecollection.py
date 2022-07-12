import numpy as np
from astropy.io import fits

from ...logger import LOGGER
from ..photometry import SED
from .directimages import DirectImages
from ..astrometry import AstroImage
from ..utilities import indices,headers

class SourceCollection(dict):
    def __init__(self,segfile,detfile,maglim=None,sedfile=None,**kwargs):
        LOGGER.info(f'Loading segmentation map {segfile}')

        self.segfile=segfile

        # load the direct images
        self.images=DirectImages.from_file(detfile,**kwargs)


        # specify a magnitude limit
        if maglim is None:
            maglim=np.inf
        self.maglim=maglim


        # open the fits files and sort out their nature
        mode='readonly'
        with fits.open(self.segfile,mode=mode) as hdus,\
             fits.open(self.obsdata.detfile,mode=mode) as hdui:

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
            if len(self) > 0:
                # load the photometry
                #self.load_photometry()

                # load the SEDs if there're present.
                if sedfile is not None:
                    self.load_seds(sedfile)
                self.sedfile=sedfile
            else:
                LOGGER.warn(f'There are no valid sources in: {self.segfile}')


    def _from_classic(self,hdus,hdui,exten=0):
        LOGGER.info(f"Loading classic segmentation image: {self.segfile}")

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

            # make a source based on the segid value
            if segid >0:
                source=SimpleSource(subimg,subseg,zeropoint=zeropoint)
            else:
                source=CompoundSource(subimg,subseg,zeropoint=zeropoint)

            # add the source if it's valid
            if source and source.mag < self.maglim:
                self[segid]=source


    def _from_mef(self,hdus,hdui):
        pass

    def unknowns(self):
        u={}
        for source in self.values():
            u.update(source.unknowns())
        return u



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
        self.obsdata.update_header(hdr)
        headers.add_stanza(hdr,'Source Properties',before='NSOURCE')


    def write_seds(self,**kwargs):

        for src in self.values():
            src.write_seds(**kwargs)

    def load_seds(self,sedfile):
        self.sedfile=sedfile
        with SEDFile(self.sedfile) as sedlib:

            for segid,source in self.items():

                for key in source.seds:
                    if isinstance(key,tuple):
                        segid,xd,yd=key
                        key=f'{segid},{xd},{yd}'
                        if key in sedlib:
                            source.set_sed(*sedlib[key],x=xd,y=yd,band=detband)
                            break
                        else:
                            LOGGER.warning(f'{key} not in SEDLIB')

                    elif isinstance(key,(int,np.integer)):
                        source.set_sed(*sedlib[key],band=detband)
                        break
                    else:
                        LOGGER.error(f'Invalid SED key type: {key} {type(key)}')
                else:
                    LOGGER.warning(f'{segid} not found in SEDLIB')

    def load_photometry(self):
        LOGGER.debug('Photometry will not load.  need to fix this for contam')
        LOGGER.info('Loading broadband photometry')

        # get some dimensions
        nbands=len(self.images)
        nsources=len(self)

        # set some variables
        lamb=np.zeros(nbands,dtype=float)
        flam=np.zeros((nbands,nsources),dtype=float)
        for i,(filename,photplam,filtfile,band) in enumerate(self.images):
            lamb[i]=photplam
            #photflam=band.photflam

            if self.mef:
                pass

            else:
                img=fits.getdata(filename,0)
                for j,(sedkey,source) in enumerate(self.items()):
                    flam[i,j]=source.instrumental_flux(img)
