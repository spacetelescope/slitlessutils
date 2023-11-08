import numpy as np

from .....logger import LOGGER
from ....utilities import headers
from .menger import menger


def optimizer(algorithm, logdamp):
    """
    Factory function to determine the optimizer

    Parameters
    ----------
    algorithm : str
        The name of the optimization

    logdamp : list, tuple, `np.ndarray`
        The meaning of this depends on the type of optimization

    Returns
    -------
    optimizer : `Optimizer`
        The optimizer
    """

    algorithm = algorithm.lower()
    if algorithm == 'golden':
        optimizer = Golden(logdamp)
    elif algorithm == 'grid':
        optimizer = GridSearch(logdamp)
    elif algorithm == 'single':
        optimizer = Single(logdamp)
    elif algorithm == 'multiple':
        optimizer = Multiple(logdamp)
    else:
        LOGGER.warn('defaulting to golden search')
        optimizer = Golden(logdamp)
    return optimizer


class Optimizer:
    """
    Base class for optimizers

    Parameters
    ----------
    logdamp : float, int, tuple, list, `np.ndarray`
       The exact nature of this variable depends on the child class

    """

    def __init__(self, logdamp):
        self.logdamp = logdamp

    def update_header(self, hdr):
        """
        Method to update a fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            The fits header to update
        """

        hdr.set('ALGORTHM', value=self.ALGORITHM,
                comment='method to optimize damping')
        hdr.set('OPTIMIZE', value=self.OPTIMIZED,
                comment='are the results optimized over the damping')
        headers.add_stanza(hdr, 'Damping Settings', before='ALGORTHM')


class GridSearch(Optimizer):
    """
    Class to implement a grid-based optimization

    Parameters
    ----------
    logdamp : list, tuple, `np.ndarray`
        A 3-element iterable, where the elements represent the starting,
        ending, and step size

    """

    OPTIMIZED = True
    ALGORITHM = 'grid'

    def __init__(self, logdamp):
        if len(logdamp) != 3:
            raise ValueError('Grid search requires a 3 element iterable')
        Optimizer.__init__(self, logdamp)

    def update_header(self, hdr):
        """
        Method to update a fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            The fits header to update
        """

        super().update_header(hdr)
        hdr.set('LOGDAMP0', value=self.logdamp[0],
                comment='lower bound of the damping')
        hdr.set('LOGDAMP1', value=self.logdamp[1],
                comment='upper bound of the damping')
        hdr.set('DLOGDAMP', value=self.logdamp[2],
                comment='threshold on the damping convergence')

    def __call__(self, matrix):
        """
        Method to call the optimizer

        Parameters
        ----------
        matrix : `su.core.modules.extract.multi.Matrix`
           The matrix operator to optimize

        Returns
        -------
        state : `su.core.modules.extract.multi.Result`
           The optimized state
        """

        # compute a grid of log damping
        logdamp = np.arange(self.logdamp[0], self.logdamp[1]+self.logdamp[2],
                            self.logdamp[2])

        maxcurv = -np.inf
        state0 = matrix.invert(logdamp[0])
        state = state1 = matrix.invert(logdamp[1])
        for i in range(1, len(logdamp)-1):
            state2 = matrix.invert(logdamp[i+1])

            curv = menger(state0.xy, state1.xy, state2.xy)
            if curv > maxcurv:
                maxcurv = curv
                state = state1

            state0 = state1
            state1 = state2

        return state


class Single(Optimizer):
    """
    Class to implement a single (unoptimized) step

    Parameters
    ----------
    logdamp : float or int
        a scalar value

    """

    OPTIMIZED = False
    ALGORITHM = 'single'

    def __init__(self, logdamp):
        if isinstance(logdamp, (tuple, list, np.ndarray)):
            raise ValueError('Single search requires a scalar')

        Optimizer.__init__(self, logdamp)

    def update_header(self, hdr):
        """
        Method to update a fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            The fits header to update
        """

        super().update_header(hdr)
        hdr.set('LOGDAMP', value=self.logdamp,
                comment='log(damping) value used')

    def __call__(self, matrix):
        """
        Method to call the optimizer

        Parameters
        ----------
        matrix : `su.core.modules.extract.multi.Matrix`
           The matrix operator to call

        Returns
        -------
        state : `su.core.modules.extract.multi.Result`
           The resultant state
        """

        state = matrix.invert(self.logdamp)

        return state


class Multiple(Optimizer):
    """
    Class to implement many (unoptimized) values

    Parameters
    ----------
    logdamp : list, tuple, `np.ndarray`
        an iterable

    """

    OPTIMIZED = False
    ALGORITHM = 'multiple'

    def __init__(self, logdamp):
        if not isinstance(logdamp, (tuple, list, np.ndarray)):
            raise ValueError('Multiple search requires an iterable')

        Optimizer.__init__(self, logdamp)

    def update_header(self, hdr):
        """
        Method to update a fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            The fits header to update
        """

        super().update_header(hdr)
        hdr.set('LOGDAMPS', value=','.join((str(ld) for ld in self.logdamp)),
                comment='log(damping) value used')

    def __call__(self, matrix):
        """
        Method to call the optimizer

        Parameters
        ----------
        matrix : `su.core.modules.extract.multi.Matrix`
           The matrix operator to iterate over

        Returns
        -------
        states : list
           The list of `su.core.modules.extract.multi.Result` objects
        """

        states = [matrix.invert(ld) for ld in self.logdamp]
        return states


class Golden(Optimizer):
    """
    Class to implement a golden-ratio based optimization

    Parameters
    ----------
    logdamp : list, tuple, `np.ndarray`
        A 3-element iterable, where the elements represent the starting,
        ending, and relative tolerance for convergence.

    Notes
    -----
    This implements the golden search method to find the point of maximum
    curvature based on the algorithm described by Culterra & Callergaro
    2020, IOP SciNotes, 1, 6.  See also:
    https://ui.adsabs.harvard.edu/abs/2020IOPSN...1b5004C/abstract

    """

    OPTIMIZED = True
    ALGORITHM = 'golden'
    PHI = (1+np.sqrt(5.))/2.

    def __init__(self, logdamp):
        if len(logdamp) != 3:
            raise ValueError('Golden search requires a 3 element iterable')
        Optimizer.__init__(self, logdamp)

    def update_header(self, hdr):
        """
        Method to update a fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            The fits header to update
        """

        super().update_header(hdr)
        hdr.set('LOGDAMP0', value=self.logdamp[0],
                comment='lower bound of the damping')
        hdr.set('LOGDAMP1', value=self.logdamp[1],
                comment='upper bound of the damping')
        hdr.set('ELOGDAMP', value=self.logdamp[2],
                comment='threshold on the damping convergence')

    def __call__(self, matrix):
        """
        Method to call the optimizer

        Parameters
        ----------
        matrix : `su.core.modules.extract.multi.Matrix`
           The matrix operator to optimize

        Returns
        -------
        state : `su.core.modules.extract.multi.Result`
           The optimized state
        """

        x1, x4, eps = self.logdamp[0], self.logdamp[1], self.logdamp[2]
        x2 = (x4+self.PHI*x1)/(1+self.PHI)
        x3 = x1+(x4-x2)
        state = [matrix.invert(x) for x in (x1, x2, x3, x4)]

        # start the Cultrera algorithm
        while (state[3].damp-state[0].damp) > eps*state[3].damp:

            # compute curvatures
            c2 = menger(state[0].xy, state[1].xy, state[2].xy)
            c3 = menger(state[1].xy, state[2].xy, state[3].xy)

            # make sure the c3 curvature is positive
            while c3 < 0:

                # swap states
                state[3] = state[2]
                state[2] = state[1]

                # select a new damping value and compute a new state
                x2 = (state[3].logdamp+self.PHI*state[0].logdamp)/(1+self.PHI)
                state[1] = matrix.invert(x2)

                # compute new curvature
                c3 = menger(state[1].xy, state[2].xy, state[3].xy)

            # now update the optimized value
            if c2 > c3:
                opt = state[1]     # optimal solution

                # swap the states
                state[3] = state[2]
                state[2] = state[1]

                # compute new damping and state
                x2 = (state[3].logdamp+self.PHI*state[0].logdamp)/(1+self.PHI)
                state[1] = matrix.invert(x2)

            else:
                opt = state[2]     # optimal solution

                # swap the states
                state[0] = state[1]
                state[1] = state[2]

                # compute new damping and state
                x3 = state[0].logdamp+(state[3].logdamp-state[1].logdamp)
                state[2] = matrix.invert(x3)

        # return the optimal value
        return opt
