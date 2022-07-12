import numpy as np

#from .polyclip import polyclip
from ..module import Module
from ...utilities import indices
from ...tables import PDTFile,PDT
from ....logger import LOGGER
from ...data.wfss import SimulatedFile


class Tabulate(Module):

    # define the pixel footprint
    DX = np.array([0,0,1,1],dtype=float)
    DY = np.array([0,1,1,0],dtype=float)


    DESCRIPTION = 'Making Tables'

    def __init__(self,ttype='pdt',nsub=5,remake=True,**kwargs):
        Module.__init__(self,self.tabulate,**kwargs)

        # set some things
        self.ttype=ttype       # which type of Table to make (PDT is normal)
        self.nsub=nsub         # the subsampling frequency (MUST BE AN INTEGER).  Typeically 4 or 5
        self.remake=remake     # boolean flag to overwrite the table

        self.pixfrac = 1.0     # the pixfrac (in a drizzle context),  DO NOT CHANGE FROM 1.0


    def __str__(self):
        s=f'Tabulation Module: \n'+super().__str__()
        return s

  
    @property
    def ttype(self):
        return self._ttype

    @ttype.setter
    def ttype(self,ttype):
        ttype=ttype.lower()
        if ttype == 'pdt':
            self._tabfunc = self.make_pdts
        elif ttype == 'rdt':
            self._tabfunc = self.make_rdts
        else:
            LOGGER.warning(f'Table type {ttype} is invalid.')
            return

        self._ttype=ttype

    @property
    def nsub(self):
        return self._nsub

    @nsub.setter
    def nsub(self,nsub):
        if not isinstance(nsub,(int,np.integer)):
            LOGGER.warning(f'Nsub should be an integer: {nsub}')
            nsub=int(nsub)

        if nsub<1:
            LOGGER.warning(f'Nsub should be >=1: {nsub}')
            nsub=1
        self._nsub=nsub

  



    def tabulate(self,data,sources,remake=None,**kwargs):
        if self._tabfunc is None:
            LOGGER.error("No function set for tabulating")
            return

        if not isinstance(remake,bool):
            remake=self.remake

        # unpack inputs
        insconf,insdata=data


        # create a wavelength grid with subsampling factor
        wav=insconf.parameters.wavelengths(nsub=self.nsub)

        with PDTFile(insdata,remake=remake,path=self.path) as h5:
            if not isinstance(h5,PDTFile):
                return h5
            outfile=h5.filename    # this will be returned by this method


            # parse the data type of the `insdata` to place content
            # in the attributes of the output file
            if isinstance(insdata,SimulatedFile):
                pass
            elif isinstance(insdata,ObservedFile):
                pass
            else:
                LOGGER.critical(f'Unknown image type: {type(insdata)}')
                outfile=None
                tabs=None

            #okay... process the file
            if outfile:
                for detname,detdata in insdata.items():
                    detconf=insconf[detname]

                    # add a detector to the file
                    h5.add_detector(detconf)

                    # process each order
                    for ordname,ordconf in detconf.items():
                        h5.add_order(ordconf)

                        # process each source
                        for segid,source in sources.items():
                            tabs=self._tabfunc(source,detdata,detconf,ordname,
                                               wav,hdf5=h5.h5order,
                                               nsub=self.nsub,**kwargs)

        # might need/want to return tables instead of filename
        return outfile



    def make_pdts(self,source,detdata,detconf,ordname,wav,hdf5=None,**kwargs):

        # HARDCODED FOR LINEARLY SAMPLED WAVELENGTHS
        
        # get the bandwidth
        wav0=np.amin(wav)
        wav1=np.amax(wav)
        nwav=len(wav)
        dwav=(wav1-wav0)/(nwav-1)

        # get the order property
        #ordconf=detconf[ordname]

        # the outputs
        pdts={}

        # the relative pixel area
        pixrat=(detdata.pixelarea/source.pixelarea)/self.pixfrac

        # get some dimensionalities
        dims=(detdata.hdr['NAXIS1'],detdata.hdr['NAXIS2'],nwav)

        # process each pixel
        for x,y,w in source:

            # coordinate in the main image
            xyd=source.image_coordinates(x,y,dtype=int)         
            
            # make an empty table
            pdt=PDT(*xyd,dims=dims,area=np.float32(source.area),
                    pixfrac=self.pixfrac,wav0=np.float32(wav0),
                    wav1=np.float32(wav1),dwav=np.float32(dwav),**kwargs)
                    
            # add some things to the PDT
            #pdt.attrs['wav0']=np.float32(wav0)
            #pdt.attrs['wav1']=np.float32(wav1)
            #pdt.attrs['dwav']=np.float32(dwav)

            # transform the pixel position and apply footprint
            xg,yg=source.xy2xy(x+self.DX,y+self.DY,detdata)

            # drizzle this pixel
            xx,yy,ll,aa=detconf.drizzle(xg,yg,wav,ordname)

            if len(xx)>0:
                # decimate this pixel.
                # in theory, this isn't needed, but if we have really small
                # direct image pixels, and/or small bandwidth (ie. large
                # NSub), then this can be important.
                aa,xx,yy,ll=indices.decimate(aa,xx,yy,ll)

                # At this point, the only effect accounted for is the
                # relative area between the grism and direct image (aa).
                # now we will include three effects:
                # 1. direct image weight (w)
                # 2. ratio of pixel areas between direct and grism (pixrat)
                # 3. wavelength sampling for integrals (dwav)
                pdt.extend(xx,yy,ll,aa*w*pixrat*dwav)


            # save this PDT
            pdts[xyd]=pdt

            # if asked to save.  save it here
            if hdf5 is not None:
                try:
                    pdt.write_hdf5(hdf5)
                except:
                    pass

        return pdts
