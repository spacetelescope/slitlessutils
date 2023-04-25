import numpy as np
import os

from ....logger import LOGGER
from ..module import Module
from ...tables import PDTFile, PDT
from ...utilities import indices, as_iterable
from ...wfss import WFSS


class Tabulate(Module):
    """
    Module to compute the dispersion tables

    inherits from `su.core.modules.Module`

    Parameters
    ----------
    orders : list or None, optional
       The spectral orders to tabulate.  If `None`, then will do all orders
       present in the configuration file.  Default is None

    ttype : str, optional
       The type of table to compute and write to disk.  Default is 'pdt'

    nsub : int, optional
       The subsampling rate.  Default is 5

    remake : bool, optional
       Flag to force remaking of the tables

    kwargs : dict, optional
       Optional keywords sent to parent class, see `su.core.modules.Module`

    Notes
    -----
    If a table exists, then the code will quietly skip over it, unless
    `remake==True`.  Otherwise, the code will make the table.

    """

    # define the pixel footprint
    DX = np.array([0, 0, 1, 1], dtype=float)
    DY = np.array([0, 1, 1, 0], dtype=float)

    DESCRIPTION = 'Tabulating'

    def __init__(self, orders=None, ttype='pdt', nsub=5, remake=True, **kwargs):
        Module.__init__(self, self.tabulate, **kwargs)

        # which type of Table to make (PDT is normal)
        self.ttype = ttype

        # the subsampling frequency (MUST BE AN INTEGER).  Typically 4 or 5
        self.nsub = nsub

        # boolean flag to overwrite the table
        self.remake = remake

        # the pixfrac (in a drizzle context)
        self.pixfrac = 1.0  # DO NOT CHANGE THIS!

        self.orders = as_iterable(orders)
        
    def __str__(self):
        s = f'Tabulation Module: \n'+super().__str__()
        return s

    @property
    def ttype(self):
        return self._ttype

    @ttype.setter
    def ttype(self, ttype):
        """
        Property setter for the table type

        Parameters
        ----------
        ttype : str
            the table type.

        Notes
        -----
        This is done as a 'property' to do some error testing.

        """

        ttype = ttype.lower()
        if ttype == 'pdt':
            self._tabfunc = self.make_pdts
        else:
            LOGGER.warning(f'Table type {ttype} is invalid.')
            return

        self._ttype = ttype

    @property
    def nsub(self):
        return self._nsub

    @nsub.setter
    def nsub(self, nsub):
        """
        Property setter for the subsampling rate

        Parameters
        ----------
        nsub : int
            the subsampling frequency

        Notes
        -----
        This is done as a 'property' to do some error testing.

        """

        if not isinstance(nsub, (int, np.integer)):
            LOGGER.warning(f'Nsub should be an integer: {nsub}')
            nsub = int(nsub)

        if nsub < 1:
            LOGGER.warning(f'Nsub should be >=1: {nsub}')
            nsub = 1
        self._nsub = nsub

    def tabulate(self, data, sources, **kwargs):
        """
        Function to perform the tabulation.

        Parameters
        ----------
        data : `WFSSCollection`
           The WFSS data to tabulate

        sources : `SourceCollection`
           The sources to tabulate

        kwargs : dict, optional
           parameters passed to the table creation.  See also
           `su.core.tables.HDF5Table`

        Returns
        -------
        result : str
           The filename of the table created

        """

        if self._tabfunc is None:
            LOGGER.error("No function set for tabulating")
            return

        # create a wavelength grid with subsampling factor
        disperser = data.disperser
        # grating=data.grating
        # grating=insconf.parameters
        # wav=insconf.parameters.wavelengths(nsub=self.nsub)

        # with PDTFile(insdata,remake=remake,path=self.path) as h5:
        with PDTFile(data, remake=self.remake, path=self.path) as h5:
            if not isinstance(h5, PDTFile):
                return h5
            outfile = h5.filename    # this will be returned by this method

            # parse the data type of the `insdata` to place content
            # in the attributes of the output file
            # if isinstance(insdata,WFSS):#SimulatedFile):
            #    pass
            # elif isinstance(insdata,WFSS): #ObservedFile):
            #    pass
            if isinstance(data, WFSS):  # SimulatedFile):
                pass
            elif isinstance(data, WFSS):  # ObservedFile):
                pass

            else:
                LOGGER.critical(f'Unknown image type: {type(insdata)}')
                outfile = None
                tabs = None

            # okay... process the file
            if outfile:

                # for detname,detdata in insdata.items():
                for detname, detdata in data.items():
                    # detconf=insconf[detname]

                    # add a detector to the file
                    # h5.add_detector(detdata.detconf)
                    h5.add_detector(detdata.config)

                    # process each order
                    orders = self.orders if self.orders else detdata.orders
                    for ordname in orders:
                        # for ordname,ordconf in detdata.config.items():
                        ordconf = detdata.config[ordname]
                        h5.add_order(ordconf)

                        # process each source
                        for segid, source in sources.items():
                            tabs = self._tabfunc(source, detdata, ordname,
                                                 disperser, hdf5=h5.h5order,
                                                 **kwargs)

        # might need/want to return tables instead of filename
        return outfile

    def make_pdts(self, source, detdata, ordname, disperser, hdf5=None, **kwargs):

        # HARDCODED FOR LINEARLY SAMPLED WAVELENGTHS

        # get the bandwidth
        wav = disperser.wavelengths(nsub=self.nsub)

        wav0 = np.amin(wav)
        wav1 = np.amax(wav)
        nwav = len(wav)
        dwav = (wav1-wav0)/(nwav-1)

        # get the order property
        # ordconf=detconf[ordname]

        # the outputs
        pdts = {}

        # the relative pixel area
        pixrat = (detdata.pixelarea/source.pixelarea)/self.pixfrac

        # get some dimensionalities
        dims = (*detdata.naxis, nwav)

        # process each pixel
        for x, y, w in source.pixels(applyltv=False, weights=True):

            # coordinate in the main image
            xd, yd = source.image_coordinates(x, y, dtype=int)

            # make an empty table
            pdt = PDT(xd, yd, dims=dims,
                      pixfrac=np.float32(self.pixfrac),
                      area=np.float32(source.area),
                      wav0=np.float32(wav0),
                      wav1=np.float32(wav1),
                      dwav=np.float32(dwav),
                      units=disperser.units, **kwargs)

            # transform the pixel position and apply footprint
            xg, yg = detdata.xy2xy(x+self.DX, y+self.DY, source.wcs, forward=False)

            # drizzle this pixel
            xx, yy, ll, aa = detdata.config.drizzle(xg, yg, ordname, wav)
            if len(xx) > 0:

                # decimate this pixel.
                # in theory, this isn't needed, but if we have really small
                # direct image pixels, and/or small bandwidth (ie. large
                # NSub), then this can be important.
                aa, xx, yy, ll = indices.decimate(aa, xx, yy, ll, dims=dims)

                # At this point, the only effect accounted for is the
                # relative area between the grism and direct image (aa).
                # now we will include three effects:
                # 1. direct image weight (w)
                # 2. ratio of pixel areas between direct and grism (pixrat)
                # 3. wavelength sampling for integrals (dwav)
                pdt.extend(xx, yy, ll, aa*w*pixrat*dwav)

                # save this PDT
                # pdts[xyd]=pdt
                pdts[(x, y)] = pdt

                # if asked to save.  save it here
                if hdf5 is not None:
                    try:
                        pdt.write_hdf5(hdf5)
                    except BaseException:
                        pass

        return pdts
