import numpy as np

from ..utilities import indices
from .hdf5table import HDF5Table


class ODT(HDF5Table):
    """
    Class for an object-dispersion table (ODT)

    inherits from `HDF5Table`

    Notes
    -----
    This is a key class in modeling the matrices

    """

    # the columns for this table
    COLUMNS = ('x', 'y', 'lam', 'val')

    def __init__(self, source, dims=None, **kwargs):
        """
        Initializer

        Parameters
        ----------
        source : `su.sources.Source`
            The source represented by this ODT

        dims : tuple or None, optional
            The dimensions of the table, passed to the `HDF5Table()`
            See that for rules of typing.  Default is None

        kwargs : dict, optional
            additional arguments passed to `HDF5Table()`

        """

        HDF5Table.__init__(self, dims=dims, **kwargs)
        self.segid = source.segid

        # collect the PDTs
        self.pdts = {}

        # collect the input pixels
        self.pixels = []

        # stuff to compute object centroid
        # self.profile = 0.0
        # self.xc = 0.0
        # self.yc = 0.0

    @property
    def name(self):
        return str(self.segid)

    def append(self, pdt):
        """
        Method to append a PDT

        Parameters
        ----------
        pdt : `su.tables.PDT`
           A pixel-dispersion table (PDT) to include in this ODT

        """

        if pdt:
            pixel = pdt.pixel
            self.pdts[pixel] = pdt
            self.pixels.append(pixel)

            # update the one-pass averages to get center
            # profile = self.profile
            # p = pdt.attrs['profile']
            # self.profile += p
            # self.xc = (self.xc*profile + pixel[0]*p)/self.profile
            # self.yc = (self.yc*profile + pixel[1]*p)/self.profile

    def decimate(self):
        """
        Method to decimate over the PDTs

        Parameters
        ----------
        None.

        Returns
        -------
        factor : float
            The decimation factor --- the fractional decrease in the size of the table
            from summing over repeated indices.  Will be `np.nan` if the table cannot
            be decimated.

        """

        factor = np.nan
        if self.pdts:

            # extract all the values, but start with the existing data
            x = self['x']
            y = self['y']
            lam = self['lam']
            val = self['val']

            for pix, pdt in self.pdts.items():
                # update
                x.extend(pdt['x'])
                y.extend(pdt['y'])
                lam.extend(pdt['lam'])
                val.extend(pdt['val'])

            if x:
                # ok... let's just save some space
                self.pdts.clear()

                # change datatypes
                x = np.asarray(x, dtype=int)
                y = np.asarray(y, dtype=int)
                lam = np.asarray(lam, dtype=int)
                val = np.asarray(val, dtype=float)

                # do the summations
                vv, xx, yy, ll = indices.decimate(val, x, y, lam, dims=self.dims)
                factor = 1 - xx.size/x.size

                # put these values in the self
                self.clear()
                self['x'] = xx
                self['y'] = yy
                self['lam'] = ll
                self['val'] = vv

        return factor
