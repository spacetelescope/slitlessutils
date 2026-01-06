from datetime import datetime

import numpy as np
from astropy.io import fits

from ....logger import LOGGER
from ...tables import PDTFile
from ...utilities import headers, indices
from ...utilities.compression import compress
from ..module import Module


class Simulate(Module):
    """
    Module to simulate WFSS images

    inherits from `su.core.modules.Module`

    Parameters
    ----------
    orders : list or None, optional
       The spectral orders to include in the image.  If `None`, then all
       of the orders in the configuration will be used.  Default is None

    addnoise : bool, optional
       Flag to add photometric noise to image.  See the dataclass:
       `su.core.wfss.config.instrumentconfig.Noise` for more information.
       Default is True

    gzip : bool, optional
       Flag to gzip the output products.  Default is True

    overwrite : bool, optional
       Flag to overwrite the output products.  Default is True

    kwargs : dict, optional
       Keywords passed to the parent class.

    """

    DESCRIPTION = 'Simulating'

    def __init__(self, orders=None, addnoise=True, gzip=True, overwrite=True,
                 **kwargs):

        Module.__init__(self, self.simulate, **kwargs)

        # save some things
        self.gzip = gzip
        self.overwrite = overwrite
        self.addnoise = addnoise
        self.orders = orders

    def __str__(self):
        s = 'Simulation Module: \n' + super().__str__()
        return s

    def simulate(self, data, sources, **kwargs):
        """
        Method to do the simulation

        Parameters
        ----------
        data : `WFSSCollection`
           The WFSS data to simulate

        sources : `SourceCollection`
           The sources to simulate

        kwargs : dict, optional
           Additional keywords used for the flat fielding

        Returns
        -------
        result : str
           The filename of the image created

        """
        if sources.sedfile is None:
            LOGGER.critical("No SEDFile found.")
            return

        # grab the time
        t0 = datetime.now()

        # open the table for reading
        with PDTFile(data, path=self.path, mode='r') as h5:

            # start a fits file
            hdul = fits.HDUList()

            # make a primary header
            phdu = fits.PrimaryHDU(header=data.config.make_header())

            # update the primary header with some useful things
            super().update_header(phdu.header)

            # put the PA_V3 in the header
            if hasattr(data, 'pa_v3'):
                phdu.header.set('PA_V3', value=data.pa_v3,
                                comment='position angle of V3-axis (deg)')
                headers.add_stanza(phdu.header, 'POINTING INFORMATION',
                                   before='PA_V3')

            # put the wavelengths in the header
            # insconf.parameters.update_header(phdu.header)
            # data.grating.update_header(phdu.header)
            # data.disperser.update_header(phdu.header)

            # put the orders in the header
            # orders=self.orders if self.orders else data.orders
            # phdu.header['NORDER']=(len(orders),'Number of orders in the image')
            # phdu.header['ORDERS']=(','.join(orders),'Order names')
            # headers.add_stanza(phdu.header,'Orders Information',
            # before='NORDER')

            # update for the reference SIAF Information
            # data.siaf.update_header(phdu.header)
            # insconf.refsiaf.update_header(phdu.header)

            # process each detector
            for detname, detdata in data.items():
                # load a detector
                # detconf=detdata.config
                # detconf=insconf[detname]
                h5.load_detector(detname)

                # create an empty Image
                dtype = detdata.extensions['science'].dtype
                sci = np.zeros(detdata.shape, dtype=dtype)

                # load a flatfield
                flatfield = detdata.config.load_flatfield(**kwargs)

                # process each order in question
                orders = self.orders if self.orders else detdata.orders
                for ordname in orders:

                    ordconf = detdata.config[ordname]

                    h5.load_order(ordname)

                    # process each source
                    for segid, source in sources.items():

                        # process each region within a source
                        for regid, region in source.items():

                            # process each pixel for each region
                            for x, y in region.pixels():
                                # NB: the previous two forloops could be
                                #     condensed into a single for-loop if
                                #     we only consider Singleton Sources
                                #     or fully decomposed sources.  Therefore,
                                #     this might be slow and there may be room
                                #     for improvement here.

                                # convert to image coordinates and load PDT
                                xd, yd = source.image_coordinates(x, y, dtype=int)

                                # convert to image coordinates and load PDT
                                pdt = h5.load_pdt(xd, yd)
                                if pdt:
                                    # extract the data
                                    xg = pdt.get('x')
                                    yg = pdt.get('y')
                                    val = pdt.get('val')
                                    wav = pdt.wavelengths()

                                    # need to apply a few things:
                                    # 1) sensitivity curve    (sens)
                                    # 2) flat field           (flat)
                                    # 3) relative pixel areas (area)
                                    # 4) source spectra       (flam)
                                    sens = ordconf.sensitivity(wav)
                                    flat = flatfield(xg, yg, wav)
                                    flam = region.sed(wav, fnu=False)

                                    # apply the corrections to the weights
                                    val *= (sens * flat * flam)

                                    # sum over wavelengths
                                    vv, yy, xx = indices.decimate(val, yg, xg,
                                                                  dims=detdata.shape)

                                    # apply pixel areas
                                    vv *= detdata.relative_pixelarea(xx, yy)

                                    # now sum into the image
                                    sci[yy, xx] += vv

                # Here we have a *NOISELESS* science image in e/s, so need to:
                # 1) create an uncertainty image
                # 2) create a data-quality image
                # 3) create/update all headers.
                # the function make_HDUs takes the noiseless sci, creates
                # all ancillary data, packages into HDUs, and adds noise
                # hdus= detdata.make_HDUs(sci,self.noisepars)

                hdus = detdata.make_HDUs(sci, addnoise=self.addnoise)

                # put the HDUs into the list, but first update some header info
                for hdu in hdus:

                    # put the SIAF info in header
                    detdata.config.siaf.update_header(hdu.header)

                    # put the flatfield info in the headers
                    if hdu.header.get('EXTNAME', '') == 'SCI':
                        flatfield.update_header(hdu.header)

                    # append the HDU
                    hdul.append(hdu)

        # compute runtime
        dt = datetime.now() - t0

        # put some times into the header
        headers.add_preamble(phdu.header,
                             runtime=(dt.total_seconds(), 'runtime in s'))
        headers.add_software_log(phdu.header)
        sources.update_header(phdu.header)    # source props

        # put the primary header in the HDUL
        hdul.insert(0, phdu)

        # write the file to disk
        filename = data.filename
        hdul.writeto(filename, overwrite=self.overwrite)

        # do we gzip the file?
        if self.gzip:
            filename = compress(filename, comptype='gz')

        return filename
