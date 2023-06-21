import h5py
import os
import numpy as np

from ...logger import LOGGER
from . import attributes
from .hdf5columns import HDF5Columns


# r	Readonly, file must exist (default)
# r+	Read/write, file must exist
# w	Create file, truncate if exists
# w- or x	Create file, fail if exists
# a	Read/write if exists, create otherwise

class HDF5File(HDF5Columns):
    """
    Class to hold an HDF5 file data

    inherits from `HDF5Columns`
    """

    # the presumed file extension
    EXT = 'h5'

    def __init__(self, wfssfile, path='tables', mode='a', remake=False):
        """
        Initializer

        Parameters
        ---------
        wfssfile : `slitlessutils.core.wfss.WFSS`
            A WFSS datatype to use

        path : str, optional
            A path to where the `HDF5File` will be written.  Default is
            'tables'

        mode : str, optional
            The file mode.  See `h5py.File()`.  Default is 'a'

        remake : bool, optional
            Flag to force to remake the table.  Default is True

        """

        self.dataset = wfssfile.dataset
        self.telescope = wfssfile.telescope
        self.instrument = wfssfile.instrument
        self.disperser = wfssfile.disperser.name
        # self.grating=wfssfile.grating.grating

        self.blocking = '' if wfssfile.blocking is None else wfssfile.blocking
        # self.pupil='' if wfssfile.pupil is None else wfssfile.pupil
        self.wfsstype = wfssfile.filetype
        self.remake = remake

        self.attrs = {}

        self.path = path
        self.mode = mode

        # if the file exists and we want to remake them, then we
        # will wipe the existing file and start over
        if self.remake and os.path.exists(self.filename):
            os.remove(self.filename)

    @property
    def filename(self):
        base = os.path.join(self.path, self.dataset)
        return f'{base}.{self.EXT}'

    def __enter__(self):
        """
        Method to implement a context manager
        """
        if self.mode == 'w' or self.mode == 'a':
            # trying to make a table
            if not os.path.exists(self.filename) or self.remake:
                # make the table
                self.h5file = h5py.File(self.filename, self.mode)

                if self.mode == 'w':
                    attributes.write(self.h5file, 'telescope', self.telescope)
                    attributes.write(self.h5file, 'instrument', self.instrument)
                    attributes.write(self.h5file, 'disperser', self.disperser)
                    # attributes.write(self.h5file,'grating',self.grating)
                    # attributes.write(self.h5file,'pupil',self.pupil)
                    attributes.write(self.h5file, 'blocking', self.blocking)
                    attributes.write(self.h5file, 'wfsstype', self.wfsstype)

            else:
                raise NotImplementedError('What is this?')

        elif self.mode == 'r':
            self.h5file = h5py.File(self.filename, self.mode)
            self.telescope = attributes.load(self.h5file, 'telescope')
            self.instrument = attributes.load(self.h5file, 'instrument')
            # self.grating=attributes.load(self.h5file,'grating')
            self.disperser = attributes.load(self.h5file, 'disperser')
            self.blocking = attributes.load(self.h5file, 'blocking')
            # self.pupil=attributes.load(self.h5file,'pupil')
            self.wfsstype = attributes.load(self.h5file, 'wfsstype')

        else:
            raise NotImplementedError(f'Invalid file mode: {self.mode}')

        return self

    def __exit__(self, et, ev, etb):
        """
        Method to close the contextmanager

        """
        self.close()

    def __del__(self):
        """
        Method to override the del operator
        """
        self.close()

    def close(self):
        """
        Method to close the file

        """

        if hasattr(self, 'h5file'):
            if self.h5file:
                if self.mode != 'r':
                    for k, v in self.attrs.items():
                        attributes.write(self.h5file, k, v)
                self.h5file.close()
        else:
            LOGGER.warning('no HDF5file is open to close')

    def __str__(self):
        return f'HDF5 File for {self.dataset}'

    def add_detector(self, detconf):
        """
        Method to add a detector to the `HDF5File`

        Parameters
        ----------
        detconf : `su.core.wfss.config.DetectorConfig`

        """

        if hasattr(self, 'h5file'):
            self.h5detector = self.h5file.require_group(detconf.name)

            if self.mode != 'r':
                attributes.write(self.h5detector, 'naxis1',
                                 self.DTYPES['x'](detconf.naxis[0]))
                attributes.write(self.h5detector, 'naxis2',
                                 self.DTYPES['y'](detconf.naxis[1]))

        else:
            LOGGER.warning("no HDF5File is loaded")

    def add_order(self, ordconf):
        """
        Method to add a spectral order to the `HDF5File`

        Parameters
        ----------
        ordconf : `su.core.wfss.config.Order`

        """

        if hasattr(self, 'h5detector'):

            self.h5order = self.h5detector.require_group(ordconf.order)

            if self.mode != 'r':
                # put some stuff here
                # attributes.write(self.h5order,'conffile',ordconf.conffile)
                attributes.write(self.h5order, 'order', ordconf.order)

        else:
            LOGGER.warning("no HDF5Detector is loaded in HDF5File")

    def load_detector(self, detname):
        """
        Method to load a detector

        Parameters
        ----------
        detname : str
            The name of the detector to load
        """
        try:
            self.h5detector = self.h5file[detname]
        except BaseException:
            LOGGER.warning(f'HDF5File does not have detector {detname}')

    def load_order(self, ordname):
        """
        Method to load a spectral order

        Parameters
        ----------
        ordname : str
            The name of the spectral order to load
        """
        try:
            self.h5order = self.h5detector[ordname]
        except BaseException:
            LOGGER.warning(f'HDF5Detector does not have order {ordname}')
