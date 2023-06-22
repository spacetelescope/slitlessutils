import numpy as np

from .hdf5table import HDF5Table

from ..utilities import indices


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

        pixel = pdt.pixel
        self.pdts[pixel] = pdt
        self.pixels.append(pixel)

    def decimate(self):
        """
        Method to decimate over the PDTs
        """

        if self.pdts:

            # extract all the values, but start with the existing data
            x = self['x']
            y = self['y']
            lam = self['lam']
            val = self['val']
            for pdt in self.pdts.values():
                x.extend(pdt['x'])
                y.extend(pdt['y'])
                lam.extend(pdt['lam'])
                val.extend(pdt['val'])

            # current size of aggregated table
            n = len(x)
            if n > 0:

                # ok... let's just save some space
                self.pdts.clear()

                # change datatypes
                x = np.array(x, dtype=int)
                y = np.array(y, dtype=int)
                lam = np.array(lam, dtype=int)
                val = np.array(val, dtype=float)

                # do the summations
                vv, xx, yy, ll = indices.decimate(val, x, y, lam, dims=self.dims)
                # m = len(xx)
                # r = float(n-m)/float(n)
                # print(f'Decimation factor: {r}')

                # put these values in the self
                self.clear()
                self['x'] = xx
                self['y'] = yy
                self['lam'] = ll
                self['val'] = vv

        else:
            pass
