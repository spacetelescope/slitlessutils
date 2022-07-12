import numpy as np
from astropy.io import fits
from datetime import datetime


from ..module import Module
from ...utilities import indices,headers,gzip
from ...tables import PDTFile,PDT
from ....logger import LOGGER
from ....info import __code__,__version__

class Simulate(Module):
    ''' WFSS Simulation Module '''

    DESCRIPTION = 'Making simulations'

    def __init__(self,addnoise=True,gzip=True,overwrite=True,**kwargs):

        Module.__init__(self,self.simulate,**kwargs)


        # save some things
        self.gzip=gzip
        self.overwrite=overwrite
        self.addnoise=addnoise


    def __str__(self):
        s=f'Simulation Module: \n'+super().__str__()
        return s

    def simulate(self,data,sources,**kwargs):
        ''' main method to simulate a WFSS File '''

        if sources.sedfile is None:
            LOGGER.critical("No SEDFile found.")
            return

        # grab the time
        t0=datetime.now()

        # grab the inputs
        insconf,insdata = data

        # open the table for reading
        with PDTFile(insdata,path=self.path,mode='r') as h5:

            # start a fits file
            hdul=fits.HDUList()

            # make a primary header
            phdu=fits.PrimaryHDU(header=insdata.phdr)

            # update the primary header with some useful things
            super().update_header(phdu.header)

            # put the wavelengths in the header
            insconf.parameters.update_header(phdu.header)

            # put the orders in the header
            orders=insconf.orders
            phdu.header['NORDER']=(len(orders),'Number of orders in the image')
            phdu.header['ORDERS']=(','.join(orders),'Order names')
            headers.add_stanza(phdu.header,'Orders Information',before='NORDER')


            # update for the reference SIAF Information
            insconf.refsiaf.update_header(phdu.header)

            # process each detector
            for detname,detdata in insdata.items():

                # load a detector
                detconf=insconf[detname]
                h5.load_detector(detname)

                # create an empty Image
                dtype=detconf.extensions['science'].dtype
                sci = np.zeros(detdata.shape,dtype=dtype)

                # load a flatfield
                flatfield = detconf.load_flatfield(**kwargs)

                # process each order in question
                for ordname,ordconf in detconf.items():
                    h5.load_order(ordname)

                    # process each source
                    for segid,source in sources.items():

                        # process each region within a source
                        for regid,region in source.items():

                            # process each pixel for each region
                            for x,y in region.pixels():
                                # NB: the previous two forloops could be
                                #     condensed into a single for-loop if
                                #     we only consider Singleton Sources
                                #     or fully decomposed sources.  Therefore,
                                #     this might be slow and there may be room
                                #     for improvement here.



                                # convert to image coordinates and load PDT
                                xd,yd=source.image_coordinates(x,y,dtype=int)
                                pdt=h5.load_pdt(xd,yd)
                                if pdt:
                                    # extract the data
                                    xg=pdt.get('x')
                                    yg=pdt.get('y')
                                    val=pdt.get('val')
                                    wav=pdt.wavelengths()

                                    # need to apply a few things:
                                    # 1) sensitivity curve    (sens)
                                    # 2) flat field           (flat)
                                    # 3) relative pixel areas (area)
                                    # 4) source spectra       (flam)
                                    sens=ordconf.sensitivity(wav)
                                    flat=flatfield(xg,yg,wav)
                                    area=detdata.pixel_area_map(xg,yg)
                                    flam=region.sed(wav,fnu=False)



                                    # apply the corrections to the weights
                                    val *= (sens*flat*area*flam)

                                    # sum over wavelengths
                                    vv,yy,xx=indices.decimate(val,yg,xg,
                                        dims=detdata.shape)

                                    # now sum into the image
                                    sci[yy,xx] += vv
                # Here we have a *NOISELESS* science image in e/s, so need to:
                # 1) create an uncertainty image
                # 2) create a data-quality image
                # 3) create/update all headers.
                # the function make_HDUs takes the noiseless sci, creates
                # all ancillary data, packages into HDUs, adds noise as requested
                #hdus= detdata.make_HDUs(sci,self.noisepars)
                hdus= detdata.make_HDUs(sci,addnoise=self.addnoise)

                # put the HDUs into the list, but first update some header info
                for hdu in hdus:

                    # put the SIAF info in header
                    detconf.siaf.update_header(hdu.header)

                    # put the flatfield info in the headers
                    if hdu.header.get('EXTNAME','')=='SCI':
                        flatfield.update_header(hdu.header)

                    # append the HDU
                    hdul.append(hdu)


        # record the end time
        t1=datetime.now()
        dt=t1-t0



        # put some times into the header
        phdu.header.set('ORIGIN',value=f'{__code__} v{__version__}',
                        after='NAXIS')
        phdu.header.set('DATE',value=t1.strftime('%Y-%m-%d'),after='ORIGIN',
                        comment='date this file was written (yyyy-mm-dd)')
        phdu.header.set('RUNTIME',value=dt.total_seconds(),after='DATE',
                        comment='run time of this file in s')

        phdu.header.set('EQUINOX',value=2000.,after='INSTRUME',
                        comment='equinox of celestial coord. system')


        # put the primary header in the HDUL
        hdul.insert(0,phdu)

        # write the file to disk
        filename=insdata.filename
        hdul.writeto(filename,overwrite=self.overwrite)

        # do we gzip the file?
        if self.gzip:
            filename=gzip.gzip(filename)

        return filename
