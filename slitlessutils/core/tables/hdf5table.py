import numpy as np
import pandas as pd
from shapely import geometry
from skimage import measure

from slitlessutils.logger import LOGGER

from ...config import Config
from . import attributes
from .hdf5columns import HDF5Columns


class HDF5Table(HDF5Columns):
    """
    Class to hold a table from the HDF5 files, such as a the pixel-
    dispersion tables (pdts) or object-dispersion tables (odts)

    inherits from `HDF5Columns`
    """

    # a lookup table for the colors of the ds9 regions
    COLORS = {'0': 'white', '+1': '#1f77b4', '+2': '#ff7f0e', '+3': '#2ca02c',
              '+4': '#d62728', '+5': '#9467bd', '+6': '#8c564b', '+7': '#e377c2',
              '+8': '#7f7f7f', '+9': '#bcbd22', '+10': '#17becf', '-1': '#aec7e8',
              '-2': '#ffbb78', '-3': '#98df8a', '-4': '#ff9896', '-5': '#c5b0d5',
              '-6': '#c49c94', '-7': '#f7b6d2', '-8': '#c7c7c7', '-9': '#dbdb8d',
              '-10': '#9edae5'}

    def __init__(self, dims=None, **kwargs):
        """
        Initializer

        Parameters
        ----------
        dims : tuple or None, optional
            The dimensionality of the file.    Default is None

        kwargs : dict, optional
            keyword/value pairs to interpreted as attributes.

        """

        self.dims = dims
        self.attrs = {k: v for k, v in kwargs.items()}

        self.compargs = Config().h5pyargs

        for column in self.COLUMNS:
            self[column] = list()

        if self.dims:
            # this presumes that the dimensions are in the same order as
            # COLUMNS.  This would break if COLUMNS variable is reordered
            for name, dim in zip(self.COLUMNS, self.dims):
                self.attrs[f'n{name}'] = self.DTYPES[name](dim)

    def __str__(self):
        return f'{self.__class__.__name__} for {self.name}'

    def __len__(self):
        return len(self[self.COLUMNS[0]])

    def __setitem__(self, k, v):
        super().__setitem__(k, list(v))

    def __mul__(self, a):
        self['val'] = self.get('val')*a
        return self

    def __rmul__(self, a):
        return self.__mul__(a)

    def __imul__(self, a):
        return self.__mul__(a)

    def __itruediv__(self, a):
        self['val'] = self.get('val')/a
        return self

    def __iter__(self):
        yield from zip(*self.values())

    @property
    def columns(self):
        """
        Property to get the column names as a tuple
        """
        return tuple(self.keys())

    @property
    def ncolumns(self):
        """
        Property to get the number of columns as an int
        """
        return len(self.keys())

    def clear(self):
        """
        Method to clear the table
        """
        for v in self.values():
            v.clear()

    def extend(self, *args):
        """
        Method to add new data to the table

        Parameters
        args : tuple
            this is all the columns, so should have the same number of
            elements as the table has columns
        """
        for column, arg in zip(self.columns, args):
            self[column].extend(arg)

    def get(self, col, dtype=None):
        """
        Method to get a column data

        Parameters
        ----------
        col : str
            The name of the column to get

        dtype : type or None
            Explicitly recast the output array.  If set to `None`, then
            use the default dtype in `HDF5Column`

        Returns
        -------
        dat : `np.ndarray`
            A numpy array of the column data

        Raises
        ------
        KeyError

        """

        if dtype is None:
            dtype = self.DTYPES[col]
        dat = np.array(self[col], dtype=dtype)
        return dat

    def wavelengths(self, lam=None):
        """
        Method to decode the wavelengths from the metadata in the
        `HDF5File` into floating-point (useful) wavelengths.

        Currently, this assumes a linear dispersion.

        Parameters
        ----------
        lam : int or None, optional
           The wavelength indices.  If None, then use the values
           in the table.

        Returns
        -------
        wav : `np.ndarray` dtype=float
           The floating-point wavelength values
        """

        if 'wav0' in self.attrs and 'dwav' in self.attrs:
            wav0 = self.attrs['wav0']
            dwav = self.attrs['dwav']

            if lam is None:
                lam = self.get('lam')

            wav = wav0+lam*dwav
            return wav

    # def compute_xyg(self):
    #    LOGGER.debug('make this use utilities')
    #
    #    xyg=np.ravel_multi_index((self.get('x'),self.get('y')),
    #                             dims=(self.attrs['nx'],self.attrs['ny']),
    #                             order='F')
    #    return xyg

    def as_pandas(self):
        """
        Method to dump this table as a pandas Dataframe

        Returns
        -------
        data : `pandas.DataFrame`
            The output data frame
        """

        data = {column: self.get(column) for column in self.COLUMNS}
        data = pd.DataFrame(data=data)
        return data

    def as_array(self):
        """
        Method to dump the table as a numpy structured array

        Returns
        -------
        data : `np.ndarray`
            A structured array of the table
        """

        # get the right dtypes
        dtype = [(column, self.DTYPES[column]) for column in self.COLUMNS]

        # create an array
        data = np.empty((len(self),), dtype=dtype)

        # fill the array
        for column in self.COLUMNS:
            data[column] = self.get(column)
        return data

    def bounding_box(self, dx=(0, 0), dy=(0, 0)):
        """
        Method to compute the bounding box from this table

        Parameters
        ----------
        dx : 2-tuple of integers, optional
            Extra padding to include in the bounding box of the form
            (left,right) padding.  Default is (0,0)

        dy : 2-tuple of integers, optional
            Extra padding to include in the bounding box of the form
            (below,above) padding.  Default is (0,0)

        Returns
        -------
        x0 : int
            lower x-value

        x1 : int
            upper x-value

        y0 : int
            lower y-value

        y1 : int
            upper y-value
        """

        if 'x' in self and 'y' in self:
            x = self.get('x')
            x0 = np.amin(x)-dx[0]
            x1 = np.amax(x)+dx[1]
            if 'nx' in self.attrs:
                x0 = np.ceil(max(x0, 0))
                x1 = np.floor(min(x1, self.attrs['nx']))
            x0, x1 = int(x0), int(x1)

            y = self.get('y')
            y0 = np.amin(y)-dy[0]
            y1 = np.amax(y)+dy[1]
            if 'ny' in self.attrs:
                y0 = np.ceil(max(y0, 0))
                y1 = np.floor(min(y1, self.attrs['ny']))
            y0, y1 = int(y0), int(y1)

        else:
            x0 = x1 = y0 = y1 = None

        return x0, x1, y0, y1

    def select(self, g):
        """
        Method to subselect particular rows of the table

        Parameters
        ----------
        g : `np.ndarray`
           The indices (or rows) to select

        Notes
        -----
        The table selection is done in place and cannot be undone without
        re-reading the table from the file.
        """

        for k in self.columns:
            self[k] = self.get(k)[g]

    def threshold(self, a):
        """
        Method to subselect rows of the table based on a threshold on
        the table's value

        Parameters
        ----------
        a : int or float
            The value to select against
        """

        g = self.get('val') >= a
        self.select(g)

    def compute_vertices(self, level=0.5, pad=3):
        """
        Method to compute the vertices that will outline the (x,y) pixels

        This uses `skimage.measure.find_contours()` on a boolean image
        (one if the pixel is inside the table, zero if it is outside).


        Parameters
        ----------
        level : int or float, optional
            The threshold to make the contours.  Default is 0.5

        pad : int, optional
            The additional padding when extracting coordinates and to
            improve the binning.  Default is 3

        Returns
        -------
        px : `np.ndarray`
            The x-coordinates of the polygon

        py : `np.ndarray`
            The y-coordinates of the polygon
        """

        if ('x' in self.COLUMNS) and ('y' in self.COLUMNS) and len(self) > 0:

            x = self.get('x')
            y = self.get('y')

            y0, y1 = np.amin(y), np.amax(y)
            x0, x1 = np.amin(x), np.amax(x)

            msk = np.zeros((y1-y0+2*pad+1, x1-x0+2*pad+1), dtype=bool)  # np.uint8)
            msk[y-y0+pad, x-x0+pad] = 1

            # contour the image
            contours = measure.find_contours(msk, level=level)

            # reset the contours
            py = contours[0][:, 0]+y0-pad
            px = contours[0][:, 1]+x0-pad

        else:
            px, py = np.array([]), np.array([])

        return px, py

    def shapelyPolygon(self, **kwargs):
        """
        Method to distill this table as a Shapely `Polygon` object

        See also `shapely.geometry.Polygon()`

        Parameters
        ----------
        kwargs : dict, optional
           Keywords to pass to `HDF5Table().compute_vertices()`

        Returns
        -------
        poly : `shapely.geometry.Polygon`
            The shapely polygon
        """

        if ('x' in self.COLUMNS) and ('y' in self.COLUMNS):
            px, py = self.compute_vertices(**kwargs)
            poly = geometry.Polygon(list(zip(px, py)))
            return poly

    def ds9region(self, order=None, mask=False, size=12, width=3, **kwargs):
        """
        Method to distill this table as a string of ds9 region information


        Parameters
        ----------
        order : str or None, optional
           The spectral order for which this is associated, which is only
           used to set the color of the ds9 region.  If None, then
           a default color of 'black' is used.  If a valid string,
           then the colors come from the class variable `COLORS`
           Default is None.

        mask : bool, optional
           A flag that the ds9 region should be interpred as a mask.
           Default is False

        size : int, optional
           The fontsize in the ds9 region.  Default is 12

        width : int, optional
           The width of the ds9 region.  Default is 3

        kwargs : dict, optional
           optional keywords sent to `compute_vertices()`

        Returns
        -------
        reg : str
           A string that contains the ds9 region

        """

        if ('x' in self.COLUMNS) and ('y' in self.COLUMNS):

            px, py = self.compute_vertices(**kwargs)

            # coord=','.join('{},{}'.format(*xy) for xy in zip(cx,cy))
            coord = ','.join(f'{x},{y}' for x, y in zip(px, py))
            font = f'helvetica {int(size)} bold'

            color = self.COLORS.get(order, 'black')

            msk = '-' if mask else ''

            reg = f'{msk}polygon({coord}) # color={color} ' +\
                f'text={{{self.name}}} edit=0 move=0 rotate=0 fixed=1 ' +\
                f'font="{font}" width={int(width)}'

            # if family not in ('helvitica','times','courier'):
            #    family='helvetica'
            # font=f'{family} {int(size)}'

            # if bold:
            #    font+=' bold'
            # if italic:
            #    font+= 'italic'

            # reg=f'{int(mask)}polygon({coord}) # color={color} '+\
            #    f'text={{{self.name}}} edit={int(edit)} move={int(move)} '+\
            #    f'rotate={int(rotate)} width={int(width)} font="{font}" '+\
            #    f'fixed={int(fixed)}'

        else:
            LOGGER.error('Cannot make region, table does not "x" and "y"')

        return reg

    def write_hdf5(self, h5):
        """
        Method to write this table to an HDF5 object

        Parameters
        ----------
        h5 : an `h5py.Group` object
           The HDF object to write to.
        """

        hd = h5.create_dataset(self.name, data=self.as_array(), **self.compargs)
        for k, v in self.attrs.items():
            attributes.write(hd, k, v)
