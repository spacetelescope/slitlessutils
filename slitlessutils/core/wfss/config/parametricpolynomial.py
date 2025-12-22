import numpy as np

from .spatialpolynomial import SpatialPolynomial


def arrayify(func):
    """
    A decorator that logs the execution time of a function.
    """
    def wrapper(self, x, y, z, **kwargs):

        x, y = np.atleast_1d(x, y)
        if np.isscalar(z) or not kwargs.get('pairwise', False):
            z = np.atleast_1d(z)
            z = z[:, np.newaxis]

        x = x.astype(float)
        y = y.astype(float)
        z = z.astype(float)

        p = np.squeeze(func(self, x, y, z, **kwargs))

        if p.ndim == 0:
            return p.item()

        return p
    return wrapper


class Polynomial(list):
    def __init__(self, data, name, maxiter=10, threshold=1e-3):
        self.name = name
        self.maxiter = maxiter
        self.threshold = threshold

        polys = []
        indices = []
        for k, v in data.items():
            tokens = k.split('_')
            if len(tokens) == 2 and tokens[0] == name:

                P = SpatialPolynomial(v)
                if P:
                    polys.append(P)
                    indices.append(int(tokens[1]))

        for i in sorted(indices):
            self.append(polys[i])

        self.order = len(self) - 1


class StandardPolynomial(Polynomial):
    def __init__(self, data, name, **kwargs):
        super().__init__(data, name, **kwargs)

        if self.order == 1:
            self.invert = self._linear
        elif self.order == 2:
            self.invert = self._quadratic
        else:
            self.invert = self._newton

    @arrayify
    def _linear(self, x, y, p):
        '''
        Analytically find t that solves:

        p = b(x, y) + m(x,y)*t
        on interval [0, 1]
        '''

        b = self[0].evaluate(x, y)
        m = self[1].evaluate(x, y)
        return np.clip((p - b) / m, 0, 1)

    @arrayify
    def _quadratic(self, x, y, p):
        '''
        Analytically find t that solves:

        p = a(x,y) + b(x,y)*t + c(x,y)*t^2

        on the interval [0,1]
        '''

        a = self[0].evaluate(x, y)
        b = self[1].evaluate(x, y)
        c = self[2].evaluate(x, y)

        const = -0.5 * b / c
        sroot = np.sqrt(const * const + (p - a) / c)

        tp = const + sroot
        tm = const - sroot

        return np.where((0 <= tp) & (tp <= 1), tp, tm)

    @arrayify
    def _newton(self, x, y, p):
        '''
        Use Newton's method with Halley update to find solution to

        p = a(x,y) + b(x,y)*t + c(x,y)*t^2 + d(x,y)*t^3 + ...

        on the interval [0,1]
        '''

        # compute polynomial coefficients
        c = np.empty((self.order + 1, x.size))
        for k, poly in enumerate(self):
            c[k, :] = poly.evaluate(x, y)

        # initialize
        t = np.full_like(p, 0.5)
        for itr in range(self.maxiter):

            # compute polynomials and derivatives
            dP2 = 0.0  # the second derivative
            dP = 0.0   # the first derivative
            P = 0.0    # the polynomial
            for i in range(self.order, -1, -1):
                dP2 = dP2 * t + 2 * dP
                dP = dP * t + P
                P = P * t + c[i, :]

            # compute a newton step
            dt = (p - P) / dP

            # update the step for a Halley tweak
            dt /= (1 + (dt / 2) * (dP2 / dP))

            # update the position
            t = t + dt

            # clip to be in range
            if np.amax(np.abs(dt)) < self.threshold:
                break

        # return and force to be in the domain
        return np.clip(t, 0, 1)

    @arrayify
    def evaluate(self, x, y, t, pairwise=False):
        '''
        Evaluate polynomial using Horner's method
        '''

        p = self[-1].evaluate(x, y)
        for k in range(self.order - 1, -1, -1):
            p = p * t + self[k].evaluate(x, y)
        return p

    @arrayify
    def deriv(self, x, y, t, pairwise=False):
        dp = 0.
        p = self[-1].evaluate(x, y)

        for k in range(self.order - 1, -1, -1):
            dp = dp * t + p
            p = p * t + self[k].evaluate(x, y)

        return dp


class LaurentPolynomial(Polynomial):
    def __init__(self, data, name, **kwargs):
        super().__init__(data, name, **kwargs)
        raise NotImplementedError("Laurent Polynomial is not finished")
