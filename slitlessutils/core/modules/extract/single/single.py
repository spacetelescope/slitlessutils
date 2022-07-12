import numpy as np
from astropy.io import fits
from datetime import datetime
#from astropy.stats import sigma_clipped_stats
#import pandas as pd

from .contam import Contamination
from .cont2 import ContamFactory
from ...module import Module
from ....tables import PDTFile
from ....utilities import indices,headers
from .....logger import LOGGER
from .....info import __code__,__version__
from .....config import Config,SUFFIXES





class Single(Module):


    DESCRIPTION = 'Single Extraction'


    def __init__(self,degree=0,mskorders=None,**kwargs):
        Module.__init__(self,self.extract,postfunc=self._combine,**kwargs)

        self.contam=Contamination(degree,mskorders,thresh=1.)



    def _combine(self,results,wfssdata,sources,**kwargs):
        LOGGER.info(f"Combining the 1d spectra")

        # output file type
        filetype='1d spectra'

        # suffix of all output data
        suffix=SUFFIXES[filetype]

        # grab the default extraction parmaeters
        defpars=wfssdata.get_parameters()

        # the below may not work if a source isn't in the first result catalog?
        #
        # aggregate the results
        result=results.pop(0)
        while results:
            r=results.pop(0)
            for segid,res in r.items():
                result[segid]['flam'].extend(res['flam'])
                result[segid]['func'].extend(res['func'])
                result[segid]['wave'].extend(res['wave'])


        # output datatypes
        dtype=[('wavelength',float),('flam',float),('func',float),('npix',int)]
        csvhdr=','.join(d[0] for d in dtype)
        csvfmt=['%f','%f','%f','%d']


        # get all the objects in this output catalog
        segids=list(result.keys())


        # make an output structure
        hdul=fits.HDUList()

        # make a primary header and fill it with useful info
        phdu=fits.PrimaryHDU()
        headers.add_preamble(phdu.header,filetype)    # basic preamble
        headers.add_software_log(phdu.header)         # about software
        wfssdata.update_header(phdu.header)           # about WFSS images
        sources.update_header(phdu.header)            # about the sources
        #defpars.update_header(phdu.header)            # default extraction pars
        super().update_header(phdu.header)            # about the CPU settings
        hdul.append(phdu)


        # get the flux funits
        CONFIG=Config()
        funits=f'{CONFIG.fluxscale} {CONFIG.fluxunits}'

        # now process each object
        for segid in segids:
            res=result.pop(segid)

            source=sources[segid]

            if hasattr(source,'parameters'):
                pars=source.parameters
            else:
                pars=defpars

            # convert to np.array to do calculation
            flam=np.array(res['flam'])
            func=np.array(res['func'])

            # find the spectral bins that these elements belong to
            ri=indices.reverse(pars.indices(np.array(res['wave'])))

            # make an output data structure
            out=np.full(len(pars),np.nan,dtype=dtype)
            out['wavelength']=pars.wavelengths()
            out['npix']=0

            # compute weighted-averages over the bins
            for ind,g in ri.items():
                wht=1./func[g]**2
                out['flam'][ind]=np.average(flam[g],weights=wht)
                out['func'][ind]=1./np.sqrt(np.sum(wht))
                out['npix'][ind]=len(g[0])



            # output the data
            hdu=fits.BinTableHDU(data=out)
            hdu.header.set('TUNIT1',value=pars.units,after='TFORM1')
            hdu.header.set('TUNIT2',value=funits,after='TFORM2')
            hdu.header.set('TUNIT3',value=funits,after='TFORM3')
            hdu.header.set('TUNIT4',value='number',after='TFORM4')
            source.update_header(hdu.header)
            pars.update_header(hdu.header)
            hdul.append(hdu)

            # output a CSV
            #np.savetxt(f'{segid}_{suffix}.csv',out,delimiter=',',fmt=csvfmt,header=csvhdr)




        # get basename of output file
        root=kwargs.get('root',__code__)
        hdul.writeto(f"{root}_{suffix}.fits",overwrite=True)

    def extract(self,data,sources,**kwargs):
        ''' extract spectra of all sources from a given WFSS file '''

        # grab the inputs
        insconf,insdata = data

        # scale of the output fluxes
        fluxscale=Config().fluxscale

        # build a contamination model for this dataset
        #cont=Contamination(insdata,sources,insconf.orders,**kwargs)
        cont=self.contam.make_model(insdata,sources)



        # create a data structure to hold the results
        results={segid:{'wave':[],'flam':[],'func':[]} for segid in sources.keys()}

        # open the PDT for reading
        with PDTFile(insdata,path=self.path,mode='r') as h5:

            for detname,detconf in insconf.items():
                h5.load_detector(detname)
                detdata=insdata[detname]

                # load the data
                sci=detdata.read_science(header=False)
                unc=detdata.read_uncertainty(header=False)
                dqa=detdata.read_dataquality(header=False)
                LOGGER.debug('get DQA in single extraction')

                # process each order
                for ordname,ordconf in detconf.items():
                    LOGGER.debug('testing extraction for +1... gotta get read orders')
                    if ordname !='+1':
                        break

                    h5.load_order(ordname)

                    # process each source
                    for segid,source in sources.items():

                        # load the object-profile table
                        opt=h5.load_opt(source)
                        xg=opt.get('x')
                        yg=opt.get('y')
                        val=opt.get('val')
                        wav=opt.get('wav')


                        #wht=cont[detname].weight(xg,yg)
                        wht=1.
                        sen=ordconf.sensitivity(wav)
                        #dq=dqa[yg,xg]
                        g=np.where((sen > 0) & (wht > 0))[0]
                        if len(g)>0:
                            xg,yg=xg[g],yg[g]
                            den=val[g]*(sen[g]*fluxscale)

                            results[segid]['flam'].extend(sci[yg,xg]/den)
                            results[segid]['func'].extend(unc[yg,xg]/den)
                            results[segid]['wave'].extend(wav[g])



        return results
