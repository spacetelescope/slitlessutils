from dataclasses import dataclass
import numpy as np
from astropy.wcs import FITSFixedWarning
import warnings
import pandas as pd
from io import StringIO
from glob import glob

from .wfss import WFSS
from ..config import InstrumentConfig
from ...utilities import headers
from ....logger import LOGGER


@dataclass
class ObservedData:
    """
    A Dataclass to configure an observed WFSS file

    Parameters
    ----------
    filename : str
       The full path to an observed file
    """

    filename: str

    def load_file(self, **kwargs):
        """
        Method to load an observed file

        Parameters
        ----------
        kwargs : dict, optional
            optional keywords passed to `su.core.wfss.data.WFSS.observed()`

        Return
        ------
        wfss : `su.core.wfss.data.WFSS`
            the WFSS data
        """

        with warnings.catch_warnings():
            warnings.simplefilter('ignore', FITSFixedWarning)
            wfss = WFSS.observed(self.filename, **kwargs)
        return wfss

    def get_parameters(self):
        """
        Method to read the dispering information

        Returns
        -------
        disperser : `su.core.wfss.config.Disperser`
            the information on the dispersive element
        """

        conf = InstrumentConfig.from_fitsfile(self.filename)
        return conf.disperser

    def npixels(self):
        """
        Method to get the number of pixels

        Returns
        -------
        npix : int
            Number of pixels in this image
        """

        conf = InstrumentConfig.from_fitsfile(self.filename)
        npix = conf.npixels()
        return npix


@dataclass
class SimulatedData:
    """
    Dataclass for Simulated WFSS data

    Parameters
    ----------
    telescope : str
        Name of the telescope

    instrument : str
        Name of the instrument

    dataset : str
        Name of the dataset

    ra : float
        RA of the field's center

    dec : float
        Dec of the field's center

    orientat : float
        orientat (or angle East of North) of the field's rotation

    disperser : str
        Name of the dispersive element

    blocking: str
        Name of the blocking filter (should be `None` for all HST instruments)

    exptime : float
        notional exposure time (in seconds).  Only used for uncertainty model.

    background : float
        notional background level (in e-/s).  Only used for uncertainty model.

    """

    telescope: str
    instrument: str
    dataset: str
    ra: float
    dec: float
    orientat: float
    disperser: str
    blocking: str
    exptime: float
    background: float

    def load_file(self, **kwargs):
        """
        Method to load a Simulated file

        Parameters
        ----------
        kwargs : dict, optional
            Optional parameters passed to `su.core.wfss.data.WFSS.observed()`

        Returns
        -------
        wfss : `su.core.wfss.data.WFSS`
            The WFSS data
        """
        wfss = WFSS.simulated(self.telescope, self.instrument, self.dataset,
                              self.ra, self.dec, self.orientat,
                              (self.disperser, self.blocking),
                              exptime=self.exptime,
                              background=self.background, **kwargs)
        return wfss

    def __post_init__(self):
        """
        Method to validate inputs

        Notes
        -----
        This should never be directly called.
        """

        if self.blocking == '':
            self.blocking = None


class WFSSCollection(dict):
    """
    An container object to keep all the WFSS data in a single place that
    emulates a `dict`.

    Notes
    -----
    This is a key object

    inherits from dict
    """

    # possible comment tokens in reading an ascii file
    COMMENTS = ('%', '#', '&', '$', '!', ';')

    def __init__(self, filename=None, filetype=None):
        """
        The initializer method

        Parameters
        ----------
        filename : str or `None`, optional
            The name of the collection

        filetype : str or `None`, optional
            The way the WFSS collection was loaded

        Notes
        -----
        This is likely not to be called by itself, but rather should
        probably be called via one of the classmethods: `from_glob()`,
        `from_list()`, `from_file()`, `from_dataframe()`, `from_wcsfile()`
        """

        self.filename = filename
        self.filetype = filetype



    def get_visits(self):
        """
        Method to get visit ids for all the files in this collection

        Returns
        -------
        visit : `np.ndarray` dtype="<U2"
            The visit identifier for each exposure
        """
        visits = [obs.load_file().visit for obs in self.values()]
        return np.asarray(visits)

    def get_orbits(self):
        """
        Method to get orbit ids for all the files in this collection

        Returns
        -------
        visit : `np.ndarray` dtype="<U2"
            The visit identifier for each exposure
        """
        orbits = [obs.load_file().orbit for obs in self.values()]
        return np.asarray(orbits)
       

    def get_pas(self,**kwargs):
        """
        Method to get the average PAs for each exposure

        Parameters
        ----------
        kwargs : dict, optional
            See `WFSSDetector.get_pa()` for the options

        Returns
        -------
        pa : `np.ndarray` dtype=float
            The position angle for each exposure
        """

        pas = [obs.load_file().get_pa(**kwargs) for obs in self.values()]
        return np.asarray(pas)
    
        
    def get_pixscls(self):
        """
        Method to get the average pixel scales for each exposure

        Returns
        -------
        px : `np.ndarray`, dtype=float
            The pixel scale in the x-axis for each exposure
        
        py : `np.ndarray`, dtype=float
            The pixel scale in the y-axis for each exposure
        """
        
        n = len(self)
        pxs = np.zeros(n, dtype=float)
        pys = np.zeros(n, dtype=float)
        for i, obs in enumerate(self.values()):
            wfss = obs.load_file()
            pxs[i], pys[i] = wfss.get_pixscl()
                
        return pxs,pys   

    def __bool__(self):
        """
        Method to overload the bool() operator

        Returns
        -------
        True or False if it has a valid telescope set
        """

        return self.telescope is not None

    def __str__(self):
        s = f'WFSS Collection: {len(self)}'
        return s

    def __setitem__(self, k, v):
        """
        Method to overload the `dict` item setter

        Parameters
        ----------
        k : str
           The name of a dataset

        v : `su.core.wfss.data.WFSS`
           A WFSS dataset.

        """

        if not isinstance(v, (SimulatedData, ObservedData)):  # WFSS):
            LOGGER.warning(f'Can only set WFSS objects')
            return

        if not hasattr(self, 'config'):

            if isinstance(v, SimulatedData):
                self.config = InstrumentConfig(v.telescope, v.instrument,
                                               (v.disperser, v.blocking))
            elif isinstance(v, ObservedData):
                self.config = InstrumentConfig.from_fitsfile(v.filename)
            else:
                msg = f"Unknown datatype: {type(v)}"
                LOGGER.warning(msg)
                raise NotImplementedError(msg)

        super().__setitem__(k, v)

    def from_index(self, i):
        """
        Method to get a dict value by index number (as opposed to name)

        Parameters
        ----------
        i : int
            Index to get

        Returns
        -------
        v : `su.core.wfss.data.WFSS`
            A WFSS data object

        Raises
        ------
        `IndexError` if index is out of bounds
        """

        keys = list(self.keys())
        return self[keys[i]]

    def get_parameters(self, index=0):
        """
        Method to get the extraction parameters

        Parameters
        ----------
        index : int, optional
            The index of the collection to use for dispersion.  Default is 0

        Returns
        -------
        disp : `su.core.wfss.config.Disperser`
            The disperser class for the parameters

        """
        wfss = self.from_index(index)
        disp = wfss.get_parameters()
        return disp

    @classmethod
    def from_glob(cls, glb):
        """
        Classmethod to load a `WFSSCollection` from a glob string

        Parameters
        ----------
        glb : str
            A string to pass to `glob.glob()`

        Returns
        -------
        obj : `WFSSCollection`
            The `WFSSCollection`

        Notes
        -----
        Only loads `ObservedData`

        """

        LOGGER.info(f'Loading from glob: {glb}')
        obj = cls(filetype='observed', filename='<glob>')
        for f in glob(glb):
            obj[f] = ObservedData(f)
        return obj

    @classmethod
    def from_list(cls, filenames, filename='<LIST>'):
        """
        Classmethod to load a `WFSSCollection` from a list

        Parameters
        ----------
        filenames : list, tuple
            A list of filenames to load

        filename : str, optional
            a name of the parent file, if the list comes from reading an
            ascii file.  Default is '<LIST>'

        Returns
        -------
        obj : `WFSSCollection`
            The `WFSSCollection`

        Notes
        -----
        Only loads `ObservedData`

        """

        LOGGER.info('Loading WFSS data from python list')

        obj = cls(filetype='observed', filename=filename)
        for f in filenames:
            obj[f] = ObservedData(f)
        return obj

    @classmethod
    def from_file(cls, obsfile):
        """
        Classmethod to load a `WFSSCollection` from an ascii file

        Parameters
        ----------
        obsfile : str
            full path to an ascii file that contains a list of the input
            images.  This file should contain one file per line and
            valid comment cards are listed above.

        Returns
        -------
        obj : `WFSSCollection`
            The `WFSSCollection`

        Notes
        -----
        Only loads `ObservedData`

        """

        LOGGER.info(f'Loading WFSS data from file list: {obsfile}')
        obj = cls(filename=obsfile, filetype='observed')
        with open(obsfile, 'r') as fp:
            for line in fp:
                line = line.strip()
                if line and line[0] not in cls.COMMENTS:
                    obj[line] = ObservedData(line)
        return obj

    @classmethod
    def from_dataframe(cls, df, **kwargs):
        """
        Classmethod to load a `WFSSCollection` from an pandas dataframe

        Parameters
        ----------
        df : `pd.DataFrame`
            A valid pandas DataFrame that has:
            1) the index is a `str` of dataset values; and
            2) the columns are given by the parameters in `SimulatedData`

        Returns
        -------
        obj : `WFSSCollection`
            The `WFSSCollection`


        Notes
        -----
        Only loads `SimulatedData`

        """

        obj = cls(filetype='simulated', **kwargs)
        for dataset, row in df.iterrows():
            obj[dataset] = SimulatedData(dataset=dataset, **row)

        return obj

    @classmethod
    def from_wcsfile(cls, wcsfile, **kwargs):
        """
        Classmethod to load a `WFSSCollection` from a comma-separated value
        (CSV) file that tabulates the simulation parameters

        Parameters
        ----------
        wcsfile : str
            The full path to an ascii file that tabulates the simulation
            parameters.  This file should have columns given by the parameters
            for a `SimulatedData` object.

        Returns
        -------
        obj : `WFSSCollection`
            The `WFSSCollection`


        Notes
        -----
        Only loads `SimulatedData`

        """

        LOGGER.info(f'Loading WFSS data from WCS csv file: {wcsfile}')

        # set some defaults (or constants)
        defaults = {'telescope': kwargs.get('telescope', 'HST'),
                    'instrument': kwargs.get('instrument', 'WFC3IR'),
                    'disperser': kwargs.get('disperser', 'G102'),
                    'blocking': kwargs.get('blocking', ''),
                    'ra': kwargs.get('ra', 42.),                     # in deg
                    'dec': kwargs.get('dec', 42.),                   # in deg
                    'orientat': kwargs.get('orientat', 0.0),         # in deg
                    'exptime': kwargs.get('exptime', 1000.),         # in sec
                    'background': kwargs.get('background', 0.0)}     # in e/s

        with open(wcsfile, 'r') as fp:
            dat = ''

            # read the CSV, but scan for a header with defaults and update
            for line in fp:
                if line[0] in cls.COMMENTS:
                    line = line[1:].strip()
                    tokens = line.split('=')
                    if len(tokens) == 2:
                        defaults[tokens[0]] = tokens[1]
                else:
                    dat += line

        # use the defaults to get the datatypes of things
        converters = {k: type(v) for k, v in defaults.items()}
        df = pd.read_csv(StringIO(dat), converters=converters, index_col='dataset')

        # update with the defaults
        for k, v in defaults.items():
            if k not in df:
                df.insert(0, k, v, allow_duplicates=False)

        # return as a dataframe
        return cls.from_dataframe(df, filename=wcsfile)

    def __iter__(self):
        """
        Method to overload the iterator

        Returns
        -------
        An iterator
        """

        for data in self.values():
            yield data.load_file()

    def items(self):
        """
        Method to overload the `dict.items()` method

        Returns
        -------
        An iterator that returns a 2-tuple of the dataset name and WFSS data,
        respectively.
        """

        for name, data in super().items():
            yield name, data.load_file()

    def npixels(self):
        """
        Method to compute total number of pixels in this `WFSSCollection`

        Returns
        -------
        npix : int
            Total number of pixels
        """
        npix = sum(v.npixels() for v in self.values())
        return npix


    def update_header(self, hdr):
        """
        Method to update a fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            A fits header to update

        Returns
        -------
        None
        """
        hdr['NWFSS'] = (len(self), 'Number of WFSS exposures')
        hdr['WFSSFILE'] = (self.filename, 'Name of WFSS file')
        headers.add_stanza(hdr, 'WFSS Observations', before='NWFSS')


if __name__ == '__main__':
    # data=WFSSCollection.from_file('files.lst')
    data = WFSSCollection.from_wcsfile('wcs.csv')

    for name, wfss in data.items():
        for detname, detdata in wfss.items():

            xg, yg = detdata.config.disperse(500, 500, '+1', np.array([10000.]))
            print(xg, yg)
