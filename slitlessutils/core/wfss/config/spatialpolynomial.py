import numba as nb
import numpy as np


@nb.njit
def horner1_scalar(c, x):
    ''' implement 1d horner method assuming x is a scalar '''
    n = len(c)
    y = c[-1]
    for k in range(n - 2, -1, -1):
        y = y * x + c[k]
    return y


@nb.njit
def horner1_vector(c, x):
    ''' implement 1d horner method assuming x is a np.ndarray '''
    n = len(c)
    y = np.full_like(x, c[-1])
    for k in range(n - 2, -1, -1):
        y = y * x + c[k]
    return y


@nb.njit
def horner2_vector(c, x, y):
    ''' implement 2d horner method assuming x, y are np.ndarray '''

    n, m = c.shape
    z = horner1_vector(c[:, -1], x)
    for k in range(m - 2, -1, -1):
        z = z * y + horner1_vector(c[:, k], x)
    return z


class SpatialPolynomial:
    def __init__(self, data):
        self.npar = len(data)

        m = int(self.triangular(self.npar))
        self.data = np.zeros((m, m))

        i = 0
        for j in range(m):
            for k in range(j + 1):
                self.data[j - k, k] = data[i]
                i += 1

        if self.npar == 1:
            self.evaluate = self._const
        else:
            self.evaluate = self._horner

    def __new__(cls, data):
        n = len(data)
        m = cls.triangular(n)
        if m.is_integer():
            return super().__new__(cls)
        else:
            return None

    @staticmethod
    def triangular(n):
        return (np.sqrt(1 + 8 * n) - 1) / 2

    def _const(self, x, y):
        return self.data[0, 0]

    def _horner(self, x, y):
        return horner2_vector(self.data, x, y)
