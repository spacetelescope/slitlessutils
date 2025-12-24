import numpy as np

from .parametricpolynomial import LaurentPolynomial, StandardPolynomial
from .sensitivity import Sensitivity


class Order:
    def __init__(self, data, order):
        # record the order
        self.order = order

        # load the trace data
        self.dispx = StandardPolynomial(data[self.order], 'DISPX')
        self.dispy = StandardPolynomial(data[self.order], 'DISPY')

        # load the wavelength data
        self.disptype = data.get('TYPE', 'grism').lower()
        if self.disptype == 'grism':
            self.displ = StandardPolynomial(data[self.order], 'DISPL')
        elif self.disptype == 'prism':
            self.displ = LaurentPolynomial(data[self.order], 'DISPL')
        else:
            msg = f'disptype is invalid ({self.disptype})'
            raise ValueError(msg)

        # load the sensitivity file
        self.load_sensitivity(data[self.order].get('SENSITIVITY', ''))

    def __str__(self):
        return f'Spectral order: {self.disptype}, {self.order}'

    @property
    def name(self):
        return str(self.order)

    @classmethod
    def from_asciifile(cls, filename, order):
        from slitlessutils.core.wfss.config import WFSSConfig
        cfg = WFSSConfig.read_asciifile(filename)
        obj = cls(cfg, order)
        return obj

    def load_sensitivity(self, sensfile, **kwargs):
        """
        Method to load the sensitivity curve from a filename.

        Parameters
        ----------
        sensfile : str
            Full path to the sensitivity file.

        kwargs : dict, optional
            keywords that get passed to the `Sensitivity` curve object.
        """
        self.sensitivity = Sensitivity(sensfile, **kwargs)

    def dispersion(self, x0, y0, wavelength=None):
        """
        Method to compute the dispersion in A/pix at some undispersed
        position and wavelength.  This is given by the derivative of the
        wavelength solution as a function of position along the spectral
        trace.

        .. math::
           \frac{d\\lambda}{dr} = \frac{d\\lambda/dt}{\\sqrt{(dx/dt)^2 + (dy/dt)^2}}
        where :math:`t` is parametric value given by:
        .. math::
           t = DISPL^{-1}(\\lambda)

        and :math:`DISPL` comes from the grismconf file.

        Parameters
        ----------
        x0 : int or float
           The undispersed x-coordinate

        y0 : int or float
           The undispersed y-coordinate

        wavelength : int, float, or None, optional
           The wavelength (in A) to compute the dispersion.  If set to
           `None`, then the bandpass-averaged wavelength is used.

        Returns
        -------
        dldr : float
           The instaneous dispersion in A/pix.

        """

        # default wavelength to the average value over the sensitivity
        if wavelength is None:
            wavelength = self.sensitivity.average_wavelength

        # compute the dispersion using:
        # t = h^-1(lambda)
        t = self.displ.invert(x0, y0, wavelength)

        # compute the derivative of the trace using:
        # r'(t) = sqrt(x'(t)^2 + y'(t)^2)
        dxdt = self.dispx.deriv(x0, y0, t)
        dydt = self.dispy.deriv(x0, y0, t)
        drdt = np.sqrt(dxdt * dxdt + dydt * dydt)

        # compute the derivative of the wavelength solution
        dldt = self.displ.deriv(x0, y0, t)

        # the ratio is dl/dr = l'(t)/r'(t)
        return dldt / drdt

    def deltas(self, x0, y0, wav):
        t = self.displ.invert(x0, y0, wav)
        dx = self.dispx.evaluate(x0, y0, t)
        dy = self.dispy.evaluate(x0, y0, t)
        return dx, dy
