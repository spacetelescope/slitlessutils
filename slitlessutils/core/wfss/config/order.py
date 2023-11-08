import numpy as np

from .parametricpolynomial import ReciprocalPolynomial, StandardPolynomial
from .sensitivity import Sensitivity
from .spatialpolynomial import SpatialPolynomial


class Order:

    def __init__(self, data, order):
        self.order = order
        self.disptype = data.get('TYPE', 'grism').lower()
        if self.disptype == 'grism':
            self.displ = StandardPolynomial()
        elif self.disptype == 'prism':
            self.displ = ReciprocalPolynomial()
        else:
            msg = f'disptype is invalid ({self.disptype})'
            raise ValueError(msg)

        self.dispx = StandardPolynomial()
        self.dispy = StandardPolynomial()

        # print("WORRY ABOUT ORDER OF COEFS")
        for k, v in data[self.order].items():
            if k.startswith('DISPX'):
                self.dispx.append(v)
            elif k.startswith('DISPY'):
                self.dispy.append(v)
            elif k.startswith('DISPL'):
                self.displ.append(v)
            elif k.startswith('TSTAR'):
                self.displ.tstar = SpatialPolynomial(v)
            elif k.startswith('SENSITIVITY'):
                self.load_sensitivity(v)
            else:
                pass

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
            wavelength = self.sens.wave_ave

        # compute the dispersion using:
        # t = h^-1(lambda)
        # x(t) = a+ b*t + c*t^2 + ...
        # x'(t) = b+ 2*c*t + ...
        # y'(t) = r+ 2*s*t + ...
        # l'(t) = u+ 2*v*t + ...
        # r'(t) = sqrt(dxdt**2 + dydt**2)
        # dl/dr = l'(t)/r'(t)

        # invert the DISPL function
        t = self.displ.invert(x0, y0, wavelength)

        # evaluate the derivatives (w.r.t. to t) of the DISPX,
        # DISPY, and DISPL functions, which came from grism conf
        dxdt = self.dispx.deriv(x0, y0, t)
        dydt = self.dispy.deriv(x0, y0, t)
        dldt = self.displ.deriv(x0, y0, t)

        # return the dispersion
        dldr = dldt / np.sqrt(dxdt * dxdt + dydt * dydt)
        return dldr

    def deltas(self, x0, y0, wav):
        """
         Method to compute the offsets  with respect to the undispersed
         position.  NOTE: the final WFSS position would be given by
         adding back the undispersed positions.

         Parameters
         ----------
         x0 : float, `np.ndarray`, or int
            The undispersed x-position.

         y0 : float, `np.ndarray`, or int
            The undispersed y-position.

         wav : `np.ndarray`
            The wavelength (in A).

         Returns
         -------
         dx : float or `np.ndarray`
            The x-coordinate along the trace with respect to the undispersed
            position.

         dy : float or `np.ndarray`
            The y-coordinate along the trace with respect to the undispersed
            position.

         Notes
         -----
         The undispersed positions (`x0` and `y0`) must be of the same
         shape, hence either both scalars or `np.ndarray`s with the same
         shape.  The dtype of the variables does not matter.  If these
         variables are arrays, then the output will be a two-dimensional
         array of shape (len(wav),len(x0)).  If they are scalars, then
         the output will be a one-dimensional array with shape (len(wav),).

         """

        if np.isscalar(x0) or np.isscalar(wav):
            t = self.displ.invert(x0, y0, wav)
            dx = self.dispx.evaluate(x0, y0, t)
            dy = self.dispy.evaluate(x0, y0, t)

        else:
            shape = (len(wav), len(x0))
            dx = np.empty(shape=shape, dtype=float)
            dy = np.empty(shape=shape, dtype=float)
            for i, (_x, _y) in enumerate(zip(x0, y0)):
                t = self.displ.invert(_x, _y, wav)
                dx[:, i] = self.dispx.evaluate(_x, _y, t)
                dy[:, i] = self.dispy.evaluate(_x, _y, t)

        return dx, dy

    @property
    def name(self):
        """
        The name of the order
        """

        return str(self.order)

    def __str__(self):
        return f'Spectral order: {self.disptype}, {self.order}'

    # @classmethod
    # def from_asciifile(cls, filename, order):
    #    from .wfssconfig import WFSSConfig
    #    data = WFSSConfig.read_asciifile(filename)
    #    obj = cls(data, order)
    #    return obj


if __name__ == '__main__':
    x = Order(1, 2)
    print(x)
