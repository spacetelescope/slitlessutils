import numpy as np

from slitlessutils.logger import LOGGER


class SpatialPolynomial(dict):
    """
    Class to implement 2d spatial polynomial of the form:

    .. math::

       p(x,y) = a + b*x + c*y + d*x^2 + e*x*y + f*y^2 + ...


    inherits from dict.  The key/value pairs are the spatial exponents
    (stored via 'Cantor pairing') and coefficients, respectively.

    """

    # used for printing
    EXP = {0: '\u2070', 1: '\u00b9', 2: '\u00b2', 3: '\u00b3', 4: '\u2074',
           5: '\u2075', 6: '\u2076', 7: '\u2077', 8: '\u2078', 9: '\u2079'}

    def __init__(self, values):
        """
        Method to initialize the object.

        Parameters
        ----------
        values : float, int, or `np.ndarray`
            The coefficients for this `SpatialPolynomial`.  Due to the
            definition of the spatial polynomial, the number of
            elements of `values` must be a triangular number.  See:

            https://en.wikipedia.org/wiki/Triangular_number

        """

        self.order = None
        if np.isscalar(values):
            self.order = 0
            self[(0, 0)] = values
        else:
            # the coefs are stored as cantor pairs.
            # https://en.wikipedia.org/wiki/Pairing_function#Inverting_the_Cantor_pairing_function]

            n = self.triangular(len(values))
            if n:
                self.order = n-1

                # old way of decoding cantor pairing
                i = 0
                for j in range(n):
                    for k in range(j+1):
                        self[(j-k, k)] = values[i]
                        i += 1
            else:
                msg = "Input must be an array whose length is a triangular number"
                LOGGER.error(msg)
                raise RuntimeError(msg)

    def __str__(self):
        s = f"Spatial polynomial of order {self.order}:\n"

        s += 'p(x,y) ='
        first = True
        for (i, j), cij in self.items():
            if cij == 0:
                continue

            sign = '+' if cij >= 0 else '-'
            dp = f'{np.abs(cij)}'
            if i == 1:
                dp += '*x'
            elif i > 1:
                dp += f'*x{self.EXP[i]}'

            if j == 1:
                dp += '*y'
            elif j > 1:
                dp += f'*y{self.EXP[j]}'

            if first and sign == '+':
                s += ' '+dp
            else:
                s += ' '+sign+' '+dp
            first = False

        # p=[f'{cij} x{self.EXP[i]} y{self.EXP[j]}' for (i,j),cij in self.items()]
        # s+='p(x,y) = '+' + '.join(p)
        return s

    @staticmethod
    def triangular(N):
        """
        static method to check if an integer is a triangular number

        Parameters
        ----------
        N : int
            Number to check

        Returns
        n : int or None
            Return the triangular length if valid or None if not
        """
        n = (np.sqrt(1+8*N)-1)/2
        if n.is_integer():
            return int(n)

    def evaluate(self, x0, y0):
        """
        Method to evaluate this spatial polynomial at a position

        Parameters
        ----------
        xy : list, tuple, `np.ndarray`
            An iterable with 2 elements, where the two elements are for
            the `x` and `y`.

        Returns
        -------
        p : float, `np.ndarray`
            The value of the polynomial
        """
        p = sum(coef*x0**i*y0**j for (i, j), coef in self.items())
        # p=sum(coef*xy[0]**i*xy[1]**j for (i,j),coef in self.items())
        return p
