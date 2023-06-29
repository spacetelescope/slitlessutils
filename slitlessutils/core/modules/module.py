import os

from ..utilities import Pool, headers
from ...logger import LOGGER


class Module:
    """
    Base class to implement the different modules --- main steps in
    spectral analysis

    Parameters
    ----------
    func : callable
        Function that will be passed to the multiprocessing module and
        will iterate over the WFSS images.  See also `su.utilities.Pool()`
        for more information.  Calling signature should look like:

            >>> def my_function(wfss,*args,**kwargs)  # doctest: +SKIP

        Where `wfss` is an instance of `WFSS()`.


    path : str, optional
        The output path for the tables, directory will be made if it does
        not already exist.  Default is 'tables'

    ncpu : int, optional
        The number of CPUs to use.  See also `su.utilities.Pool()` for
        more information. Default is `None`

    postfunc : callable, optional
        Function that will be called on the results of the `func()` call.
        See also `su.utilities.Pool()` for more information.  Calling
        signature should look like:

            >>> def my_function(results,data,*args,**kwargs)  # doctest: +SKIP

        where `results` are the results from the `func()` call. Default
        is `None

    multiprocess : bool, optional
        Boolean flag to use the multiprocessor.  Default is True.


    Notes
    -----
    It is unlikely this class is directly instantiated, but rather it is
    primarily used to subclass other modules.

    """

    def __init__(self, func, path='tables', ncpu=None, postfunc=None,
                 multiprocess=True, **kwargs):

        self.multiprocess = multiprocess
        self.ncpu = ncpu
        self.func = func
        self.postfunc = postfunc
        self.path = path
        if not os.path.isdir(self.path):
            os.mkdir(self.path)

    def __str__(self):
        return '\n'.join([f'  path: {self.path}', f'  ncpu: {self.ncpu}'])

    def __call__(self, data, sources, **kwargs):
        """
        Method to call this module

        Parameters
        ----------
        data : `WFSSCollection`
           The WFSS data

        source : `SourceCollection`
           The sources

        kwargs : dict, optional
           Optional data passed to the callable functions

        Returns
        -------
        out : list
           The results for each datum in the `WFSSCollection`.  The are
           certain to be in order.

        Notes
        -----
        This is the primary entry point.

        """

        # a description of this process
        if hasattr(self, 'DESCRIPTION'):
            desc = self.DESCRIPTION
        else:
            desc = ''

        if len(data) == 0:
            LOGGER.warning(f'There are no WFSS images to process {desc}')
            return

        if len(sources) == 0:
            LOGGER.warning(f"There are no sources to process {desc}")
            return

        if self.multiprocess:

            # build a custom pool object:
            pool = Pool(self.func, ncpu=self.ncpu, desc=desc)

            # process each wfssdatum
            results = pool(data, sources, total=len(data), **kwargs)

            # if necessary postprocess
            if callable(self.postfunc):
                results = self.postfunc(results, data, sources, **kwargs)
        else:
            results = self.func(data, sources, **kwargs)

        return results

    @property
    def path(self):
        return self._path

    @path.setter
    def path(self, path):
        try:
            if not os.path.isdir(path):
                os.mkdir(path)

            self._path = path
        except BaseException:
            LOGGER.warn(f'Cannot make directory: {path}')

    def update_header(self, hdr):
        """
        Method to update a fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header()`
            The fits header to update
        """

        ncpu = self.ncpu if self.ncpu else 1
        hdr['NCPU'] = (ncpu, 'Number of CPUs used')
        hdr['HDF5PATH'] = (self.path, 'path to HDF5 tables')
        headers.add_stanza(hdr, 'Module Settings', before='NCPU')

    @staticmethod
    def as_set(data):
        """
        Staticmethod to convert an input into a set

        Parameters
        ----------
        data : any type
            the data to 'setify'.

        Notes
        -----
        This can take any type, as it carefully deals with strings
        """

        if isinstance(data, str):
            out = set()
            out.add(data)
        else:
            try:
                out = set(data)
            except BaseException:
                out = set((data,))
        return out
