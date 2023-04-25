import numpy as np
from .spatialpolynomial import SpatialPolynomial


class ParametricPolynomial(list):
    """
    Class to implement a nested polynomial

    Parameters
    ----------
    maxiter : int, optional
        Maximum number of iterations for the Newton-Raphson Method.  default is 10

    threshold : float, optional
        Threshold for convergence of the Newton-Raphson Method.  default is 1e-3
    """

    def __init__(self, maxiter=10, threshold=1e-3):

        # set a default for the parametric order
        self.order = 0
        self.invert = lambda x, y, f: None    # a default functional inversion

        # set some defaults
        self.maxiter = maxiter
        self.threshold = threshold

    def append(self, coefs):
        """
        Method to include another parametric order (overrides the
        list.append method).

        Parameters
        ----------
        coefs : list, tuple, or `np.ndarray`
            The coefficients to pass to the `SpatialPolynomial` object.
            The length of this iterable must be a triangular number.

        """

        poly = SpatialPolynomial(coefs)
        if poly:
            super().append(poly)
            self.order += 1

            if self.order == 1:
                self.invert = self._first
            else:
                self.invert = self._nth

    def coefs(self, x, y):
        """
        Method to evaluate all the `SpatialPolynomial` coefficients
        for a given position.

        Parameters
        ----------
        x : int or float
            The x-coordinates for the `SpatialPolynomial`s

        y : int or float
            The y-coordinates for the `SpatialPolynomial`s

        Returns
        -------
        coefs : list
            The coefficients for each parameteric order

        """

        coefs = [poly.evaluate(x, y) for poly in self]
        return coefs


class StandardPolynomial(ParametricPolynomial):
    """
    Class to implement a nested parametric polynomial of the form:

    .. math::

       p(t|x,y) = a(x,y) + b(x,y)*t + c(x,y)*t^2 + ....

    where :math:`a(x,y)`, :math:`b(x,y)`, and so on are `SpatialPolynomial`s.

    inherits from `ParametricPolynomial`

    """

    def __init__(self, **kwargs):
        ParametricPolynomial.__init__(self, **kwargs)

    def evaluate(self, x, y, t):
        """
        Method to evaluate the polynomial at some position and parameter

        Parameters
        ----------
        x : float or int
            The x-spatial coordinate to pass to the `SpatialPolynomial`s

        y : float or int
            The y-spatial coordinate to pass to the `SpatialPolynomial`s

        t : float, int, or `np.ndarray`
            The parameter to evaluate this `ParametricPolynomial`

        Returns
        -------
        f : float, or `np.ndarray`
            The value of the polynomial
        """

        return sum(p.evaluate(x, y)*t**i for i, p in enumerate(self))

    def deriv(self, x, y, t):
        """
        Method to evaluate the derivative of the polynomial at some
        position and parameter

        Parameters
        ----------
        x : float or int
            The x coordinates for the `SpatialPolynomial`s

        y : float or int
            The y coordinates for the `SpatialPolynomial`s

        t : float, int, or `np.ndarray`
            The parameter to evaluate this `ParametricPolynomial`

        Returns
        -------
        dfdt : float, or `np.ndarray`
            The value of the derivative of the polynomial
        """

        return sum(p.evaluate(x, y)*i*t**(i-1) for i, p in enumerate(self[1:], start=1))

    def _first(self, x, y, f):
        """
        Helper method to analytically invert a 1st order polynomial

        Parameters
        ----------
        x : float or int
           The x-coordinates for the `SpatialPolynomial`s

        y : float or int
           The y-coordinates for the `SpatialPolynomial`s

        f : float
           The value of the polynomial

        Returns
        -------
        t : float
           The parameter that gives the value `f`.
        """

        coefs = self.coefs(x, y)
        return (f-coefs[0])/coefs[1]

    def _nth(self, x, y, f):
        """
        Helper method to implement arbitrary order root-finding via the
        Newton-Raphson algorithm with a Halley correction.

        Parameters
        ----------
        x : float or int
           The x-coordinates for the `SpatialPolynomial`s

        y : float or int
           The y-coordinates for the `SpatialPolynomial`s

        f : float
           The value of the polynomial

        Returns
        -------
        t : float
           The parameter that gives the value `f`

        Notes
        -----
        For a second-order polynomial, this will terminate in one step and will
        be exactly correct.
        """

        p = np.polynomial.Polynomial(self.coefs(x, y))
        dpdt = p.deriv()
        d2pdt2 = dpdt.deriv()

        # compute the polynomial coefficients
        # coefs=np.asarray(self.coefs(x,y))[::-1]

        # compute the derivative coeffs
        # i=np.arange(-self.order,0)                    # for recip
        # dcoefs=coefs[:-1]*i
        # d2coefs=dcoefs*(i-1)

        # initialize and start the Newton-Raphson method
        convmsg = 'reach maximum number of iterations'
        t = np.full_like(f, 0.5, dtype=float)
        for itn in range(self.maxiter):

            # prepare for a NR step
            funct = p(t)
            deriv = dpdt(t)

            # A Newton-Raphson step
            dt = (funct-f)/deriv

            # A Halley update
            second = d2pdt2(t)
            dt /= (1-0.5*dt*(second/deriv))

            # apply the step and force range
            t = np.clip(t-dt, 0., 1.)

            # check for early convergence
            if np.amax(np.abs(dt)) < self.threshold:
                convmsg = 'reached required tolerance'
                break

        return t


class ReciprocalPolynomial(ParametricPolynomial):
    """
    Class to implement a nested parametric polynomial of the form:

    .. math::

       p(t|x,y) = a(x,y) + b(x,y)/(t-t*) + c(x,y)/(t-t*)^2 + ....

    where :math:`a(x,y)`, :math:`b(x,y)`, and so on are `SpatialPolynomial`s.

    inherits from `ParametricPolynomial`
    """

    def __init__(self, **kwargs):
        ParametricPolynomial.__init__(self, **kwargs)
        self.tstar = SpatialPolynomial([0.])

    def evaluate(self, x, y, t):
        """
        Method to evaluate the polynomial at some position and parameter

        Parameters
        ----------
        x : float or int
            The x-spatial coordinate to pass to the `SpatialPolynomial`s

        y : float or int
            The y-spatial coordinate to pass to the `SpatialPolynomial`s

        t : float, int, or `np.ndarray`
            The parameter to evaluate this `ReciprocalPolynomial`

        Returns
        -------
        f : float, or `np.ndarray`
            The value of the polynomial
        """

        tstar = self.tstar.evaluate(x, y)
        omega = 1./(t-tstar)
        return sum(p.evaluate(x, y)*omega**i for i, p in enumerate(self))

    def deriv(self, x, y, t):
        """
        Method to evaluate the derivative of the polynomial at some
        position and parameter

        Parameters
        ----------
        x : float or int
            The x coordinates for the `SpatialPolynomial`s

        y : float or int
            The y coordinates for the `SpatialPolynomial`s

        t : float, int, or `np.ndarray`
            The parameter to evaluate this `ParametricPolynomial`

        Returns
        -------
        dfdt : float, or `np.ndarray`
            The value of the derivative of the polynomial
        """

        tstar = self.tstar.evaluate(x, y)
        omega = 1./(t-tstar)
        return sum(p.evaluate(x, y)*i*omega**(i-1) for i, p in enumerate(self[1:], start=1))

    def _first(self, x, y, f):
        """
        Helper method to analytically invert a 1st order polynomial

        Parameters
        ----------
        x : float or int
           The x-coordinates for the `SpatialPolynomial`s

        y : float or int
           The y-coordinates for the `SpatialPolynomial`s

        f : float
           The value of the polynomial

        Returns
        -------
        t : float
           The parameter that gives the value `f`.
        """

        coefs = self.coefs(x, y)

        t = coefs[1]/(f-coefs[0])+self.tstar.evaluate(x, y)
        return t

    def _nth(self, x, y, f):
        """
        Helper method to implement arbitrary order root-finding via the
        Newton-Raphson algorithm with a Halley correction.

        Parameters
        ----------
        x : float or int
           The x-coordinates for the `SpatialPolynomial`s

        y : float or int
           The y-coordinates for the `SpatialPolynomial`s

        f : float
           The value of the polynomial

        Returns
        -------
        t : float
           The parameter that gives the value `f`.
        """

        # compute the polynomial coefficients
        coefs = np.asarray(self.coefs(x, y))[::-1]

        # compute the derivative coeffs
        i = np.arange(1-self.order, 0)                    # for recip
        dcoefs = coefs[:-1]*i
        d2coefs = dcoefs*(i-1)

        # p=lambda t: np.polyval(coefs,1./(t-tstar))
        # dpdt=lambda t: np.polyval(dcoefs,1./(t-star))/(t-tstar)**2
        # d2pdt2=lambda t: np.polyval(d2coefs,1./(t-star))/(t-tstar)**3

        # compute the bias
        tstar = self.tstar.evaluate(x, y)                # for recip

        # initialize and start the Newton-Raphson method
        convmsg = 'reach maximum number of iterations'
        t = np.full_like(f, 0.5, dtype=float)
        for itn in range(self.maxiter):

            # prepare for a NR step
            omega = 1./(t-tstar)                        # for recip
            funct = np.polyval(coefs, omega)             # for recip
            deriv = np.polyval(dcoefs, omega)*omega**2   # for recip

            # A Newton-Raphson step
            dt = (funct-f)/deriv

            # A Halley update
            second = np.polyval(d2coefs, omega)*omega**3
            dt /= (1-0.5*dt*(second/deriv))

            # apply the step and force range
            t = np.clip(t-dt, 0., 1.)

            # check for early convergence
            if np.amax(np.abs(dt)) < self.threshold:
                convmsg = 'reached required tolerance'
                break

        return t
