import os

import numpy as np

from ....logger import LOGGER
from .flatfield import load_flatfield
from .order import Order
from .pom import load_pom


class WFSSConfig(dict):

    COMMENTS = ("%", '#', '!', ';')
    KEYWORDS = ('DISPX', 'DISPY', 'DISPL', 'TSTAR', 'SENSITIVITY')

    def __init__(self, conffile, fwcpos=None, band=None, orders=None):
        """

        Parameters
        ----------
        conffile : str
            Filename for the configuration file.  This should be in
            the `grismconf` format.

        fwcpos : float or None, optional
            The 'filter wheel position', which is an angle (in degree)
            of the filter wheel.  Used only for JWST/NIRISS.  A value of
            'None' is interpreted as 0 deg.  Default value is None.

        band : str or None, optional
            The name of the direct image filter, which may dictate any
            wedge offsets.  Used only for WFC3/IR.  A value of 'None' is
            interpreted as blank.  Default value is None.

        """

        # set somethings about config file
        self.confpath = os.path.dirname(conffile)
        self.conffile = conffile

        # set some defaults
        self.ffname = None
        self.set_rotation(0.)
        self.set_shift(0., 0.)

        # read the data
        self.data = self.read_asciifile(self.conffile)

        # now parse the data

        # if there is an angle for the Filter Wheel
        if self.data.get("FWCPOS_REF"):
            self.set_fwcpos(fwcpos)

        # check for wedges in the filters
        if self.data.get(f'WEDGE_{band}'):
            self.set_shift(*self.data[f'WEDGE_{band}'])

        # look for a flatfield
        if self.data.get('FFNAME'):
            self.ffname = os.path.join(self.confpath, self.data['FFNAME'])

        # get pick-off mirror (POM) size
        self.pom = load_pom(**self.data)

        # load the orders
        for order, data in self.data.items():
            if isinstance(data, dict):
                if (isinstance(orders, (tuple, list)) and order in orders) or \
                   (not orders):

                    self[order] = Order(self.data, order)

        # some error checking
        if len(self) == 0:
            LOGGER.warning('No spectral orders are loaded.')

    def set_fwcpos(self, fwcpos):
        """
        Method to set a filter wheel rotation

        Parameters
        ----------
        fwcpos : float
            Filter Wheel position in deg

        Returns
        -------
        None
        """
        self.set_rotation(fwcpos-self.data['FWCPOS_REF'])

    def set_rotation(self, theta):
        """
        Method to set a rotation in the spectral trace.

        Parameters
        ----------
        theta : float
            Rotation angle in degrees

        Returns
        -------
        None
        """
        self.theta = np.radians(theta)
        self.cos = np.cos(theta)
        self.sin = np.sin(theta)

    def add_rotation(self, dtheta):
        """
        Method to update the current spectral dispersion  rotation angle

        Parameters
        ----------
        theta : float
            Rotation angle in degrees

        Returns
        -------
        None
        """
        self.set_rotation(self.theta+np.radians(dtheta))

    def set_shift(self, dx, dy):
        """
        Method to set a shift in the spectral trace.

        Parameters
        ----------
        dx : float
            Shift in x in pixels.

        dy : float
            Shift in y in pixels.

        Returns
        -------
        None
        """
        self.xshift = dx
        self.yshift = dy

    def add_shift(self, dx, dy):
        """
        Method to update the shift in the spectral trace.

        Parameters
        ----------
        dx : float
            Shift in x in pixels.

        dy : float
            Shift in y in pixels.

        Returns
        -------
        None
        """

        self.xshift += dx
        self.yshift += dy

    def load_flatfield(self, unity=False):
        """
        Method to read a flat field from the configuration file.

        Parameters
        ----------
        unity : boolean, optional
            Flag to return a unity flat.  Default is False

        Returns
        -------
        flatfield : `np.ndarray`
            The two-dimensional flat field
        """

        if unity or (self.ffname is None):
            ff = load_flatfield()
        else:
            ff = load_flatfield(self.ffname)
        return ff

    def dispersion(self, x0, y0, order, wavelength=None):
        """
        Method to determine the dispersion (ie. the instaneous derivative
        of the wavelength as a function of position along the trace) at
        a given position and spectral order.

        See also `slitlessutils.core.wfss.config.order.py`.

        Parameters
        ----------
        x0 : float
            The x-position for which to compute the dispersion.

        y0 : float
            The y-position for which to compute the dispersion.

        order : str
            The name of the spectral order

        wavelength : float or None, optional
            The wavelength to compute the dispersion.  If set as 'None',
            then interpreted as the 'PHOTPLAM' for the spectroscopic
            sensitivity.  Default is None.

        Returns
        -------
        dispersion : float
            The dispersion in A/pix.

        """
        return self[order].dispersion(x0, y0, wavelength=wavelength)

    def disperse(self, x0, y0, order, wav):
        """
        Method to convert an undispersed position into a WFSS image
        position using the calibrated trace and dispersion polynomials.

        Parameters
        ----------
        x0 : float or `np.ndarray`
            The undispersed x-position (in pixels).

        y0 : float or `np.ndarray`
            The undispersed y-position (in pixels).

        order : str
            The name of the spectral order

        wav : float or `np.ndarray`
            The wavelengths (in A) to disperse.

        Returns
        -------
        x : float or `np.ndarray`
            The WFSS image pixel coordinate.

        y : float or `np.ndarray`
            The WFSS image pixel coordinate.

        Notes
        -----
        The coordinates `x0` and `y0` must be of the same datatype. For
        example, if one is a `np.ndarray`, then they both must be
        `np.ndarray` (and of the same shape, but dtype is not relevant).
        """

        dx, dy = self[order].deltas(x0, y0, wav)

        x = x0+self.xshift+self.cos*dx+self.sin*dy
        y = y0+self.yshift-self.sin*dx+self.cos*dy

        return x, y

    def sensitivity(self, order, wav, **kwargs):
        """
        Method to compute the sensitivity for a given order at a
        selection of wavelengths.

        Parameters
        ----------
        order : str
            The spectral order name

        wav : float or `np.ndarray`
            The wavelengths for which the sensitivity is computed

        kwargs : dict, optional
            Keywords passed to `slitlessutils.core.wfss.config.order.py`,
            which largely govern the interpolation.

        Returns
        -------
        sensitivity : float or `np.ndarray`
            The sensitivity (in units of the sensitivity file).  Will be of
            same datatype and/or shape as `wav`.

        """

        return self[order].sensitivity(wav, **kwargs)

    @staticmethod
    def read_asciifile(conffile):
        """
        Method to read grism configuration file in the ascii format into a
        dict with some keywords specially parsed

        Parameters
        ---------
        conffile : str
            Full path to a configuration file

        Returns
        -------
        data : dict
            A dictionary of keyword/value pairs from the file

        """

        # main output structure:
        data = {}

        # open the file
        with open(conffile) as fp:

            # parse the file, line-by-line and look for commented lines
            for line in fp:
                line = line.strip()
                tokens = line.split()
                if tokens and tokens[0][0] not in WFSSConfig.COMMENTS:
                    key = tokens.pop(0)

                    # here 'key' will become the dict key and is the
                    # left-most column in the config file.  'tokens' will
                    # be further parsed to determine what to do
                    n = len(tokens)
                    if n == 0:
                        # if no tokens, then this is a "order marker" and an
                        # order will be interpreted as a sub-dict
                        order = key.split('_')[1]
                        data[order] = {}
                    else:

                        # something to hold the parsed values
                        values = []
                        for token in tokens:
                            # try casting as int, then float, then str
                            try:
                                value = int(token)
                            except ValueError:
                                try:
                                    value = float(token)
                                except ValueError:
                                    value = token

                            # save the recast value
                            values.append(value)
                        if len(values) == 1:

                            # if there is only one value, the pop it. else
                            # convert to a np array
                            if isinstance(values[0], str):
                                values = values[0]
                            else:
                                values = np.asarray(values)
                        else:
                            values = np.asarray(values)

                        # now further parse the keyword, to remove the
                        # order name
                        k = key.split('_')
                        if k[0] in WFSSConfig.KEYWORDS and k[1] in data:
                            if key.startswith('SENSITIVITY'):
                                # if the keyword indicates the sensitivity
                                # curve, then check the file exists
                                confpath = os.path.dirname(conffile)
                                newkey = 'SENSITIVITY'
                                values = os.path.join(confpath, values)
                            else:
                                order = k[1]
                                k.remove(order)
                                newkey = '_'.join(k)

                            # record the key/value pair for this order
                            data[order][newkey] = values
                        elif key.startswith('WEDGE') and len(values) == 2:
                            # check for a wedge-offset table
                            data[key] = tuple(values)
                        else:
                            # fall back condition
                            data[key] = values

        return data
