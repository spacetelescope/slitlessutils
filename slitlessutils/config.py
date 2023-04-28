
import configparser
import glob
import json
import os

from .logger import LOGGER
from .core.utilities import headers

from os import environ
from pathlib import Path
from tarfile import open as tar_open
from tempfile import gettempdir
from astropy.utils.data import download_file


def download_config_files(extract_path=gettempdir()):
    '''Method to automatically download the config files.'''

    # Cast to pathlib path object (in case we got a string)
    extract_path = Path(extract_path)

    # Download the tar file
    config_url = 'https://stsci.box.com/shared/static/1l2ploh74s4viltu491elxh82q6ms96w.gz'
    config_archive = Path(download_file(config_url, cache=True))

    with tar_open(config_archive) as tar:
        # Only extract the files that are missing
        for file in tar:
            if not (extract_path / Path(file.name)).exists():
                tar.extract(file, path=extract_path)

    # Prepare environment required for slitlessutils to be imported
    environ['slitlessutils_config'] = str(extract_path / 'slitlessutils_config')


# file suffixes.  Probably shouldn't ever change these, but here they are:
SUFFIXES = {'1d spectra': 'x1d',
            '2d spectra': 'x2d',
            '3d spectra': 'x3d',
            'L-curve': 'lcv',
            'group': 'grp',
            'matrix': 'mat',
            'wfssimage': 'flt'}


class Parameter:
    """
    A class to contain a configuration parameter


    Parameters
    ----------
    keyword : str
        The fits-header based keyword name (ie. must be <=8 char long)

    value : any type (must implement str())
        The keyword value

    comment : str
        The fits header comment card

    editable : bool, optional
        Flag that this parameter can be editted

    """

    def __init__(self, keyword, value, comment, editable=True):
        self.keyword = keyword
        self.value = value
        self.comment = comment
        self.editable = editable

    @classmethod
    def from_config(cls, p):
        """
        Classmethod to load config parameter from a config file

        Parameters
        ----------
        p : `configparser.Section`
            The section from the config file

        Returns
        -------
        obj : `Parameter`
            The output parameter
        """

        try:
            value = int(p['value'])
        except BaseException:
            try:
                value = float(p['value'])
            except BaseException:
                value = p['value']

        obj = cls(p['keyword'], value, p['comment'])

        return obj

    def __iter__(self):
        """
        Method to implement an iter (so can be cast as a `dict`)
        """
        yield ('value', self.value)
        yield ('keyword', self.keyword)
        yield ('comment', self.comment)

    def update_header(self, hdr):
        """
        Method to update an `astropy.io.fits.Header`

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            The fits header to update

        """
        hdr.set(self.keyword, value=str(self.value), comment=self.comment)

    def __str__(self):
        return f'{self.keyword:8} {self.value}'


class Config(dict):
    """
    A singleton class to hold and establish the global configuration

    Parameters
    ----------
    kwargs : dict, optional
       A dictinary of keyword/value pairs to override default values.

    """

    # a list of required keywords:
    REQUIRED = ('fluxscale', 'fluxunits', 'compression', 'compression_opts')

    # the default config filename:
    DEFAULTS = 'defaults.cfg'

    # name of the environment variable
    ENVVAR = 'slitlessutils_config'

    # enable the singleton
    _instance = None

    def __new__(cls, conffile=None, **kwargs):
        if not cls._instance:

            # make a new config object
            cls._instance = super().__new__(cls)

            # if no conffile is set, then grab the default
            if not isinstance(conffile, str):
                conffile = os.path.join(os.environ[cls.ENVVAR], cls.DEFAULTS)

            # load the config
            config = configparser.ConfigParser()
            config.read(conffile)
            cls._instance.conffile = conffile

            # load each section (which is a different parameter)
            for k in config.sections():
                cls._instance[k] = Parameter.from_config(config[k])

            # if confpath isn't specified, then grab it from the os
            if 'confpath' not in cls._instance:
                cls._instance.load_confpath(os.environ[cls.ENVVAR])

            # check that the required elements are there
            for req in cls._instance.REQUIRED:
                if req not in cls._instance:
                    LOGGER.warning(f"Missing keyword: {req}")
                    break
            else:
                # ok, this is a valid object
                pass

        # update any keywords
        for k, v in kwargs.items():
            if k in cls._instance:
                cls._instance[k].value = v

        return cls._instance

    def load_confpath(self, confpath):
        """
        Method to load a configuration file path

        Parameters
        ----------
        confpath : str
            The path of the configuration files.

        Notes
        -----
            Will require a `version.json` to exist in that directory.
        """

        versfile = os.path.join(confpath, 'version.json')

        if os.path.exists(versfile):

            self['confpath'] = Parameter('confpath', confpath, '')
            self.versfile = versfile
            with open(self.versfile, 'r') as fp:
                data = json.load(fp)
                self['confvers'] = Parameter('confvers', data[0]['version'],
                                             'config version', editable=False)
                self['confdate'] = Parameter('confdate', data[0]['date'],
                                             'config date', editable=False)
        else:
            LOGGER.error(f'Invalid confpath, missing the version file.')

    def __getattr__(self, k):
        """
        Method to override the getattr to enable attribute-like access

        Parameters
        ----------
        k : str
            Name of the attribute

        Returns
        -------
        v : any type
            The value of the attribute
        """
        if k in self:
            return self[k].value

    def __setattr__(self, k, v):
        """
        Method to override the setattr to enable attribute-like access

        Parameters
        ----------
        k : str
            Name of the attribute

        v : any type
            The value of the attribute
        """

        if k in self:
            if self[k].editable:
                if k == 'confpath':
                    self.log_confpath(v)
                else:
                    self[k].value = v
            else:
                LOGGER.error(f"Cannot set {k}={v}")
        else:
            super().__setattr__(k, v)

    def update_header(self, hdr):
        """
        Method to update a fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            The fits header to update
        """

        for k, v in self.items():
            v.update_header(hdr)

    def write(self, filename):
        """
        Method to write config file to disk as a *.cfg type

        Parameters
        ----------
        filename : str
            The name of the config file

        Notes
        -----
        will check the extension is one of (cfg, cnf, conf, config)

        """

        ext = os.path.splitext(filename)[-1][1:]
        if ext not in ('cfg', 'cnf', 'conf', 'config'):
            LOGGER.warning(f'Invalid extension: {ext}')
            return

        config = configparser.ConfigParser()
        for k, v in self.items():
            if v.editable:
                config[k] = dict(v)
        with open(filename, 'w') as f:
            config.write(f)

    def get_reffile(self, conffile, path='', **kwargs):
        """
        Method to get a refernece file from the local directory tree

        Parameters
        ----------
        conffile : str
            The filename of the configuration file

        path : str, optional
            The subpath within the configuration tree to find the file.
            default is ''

        kwargs : dict, optional
            extra parameters passed to `glob.glob()`

        Returns
        -------
        filename : str
            The full path to the configuration file.  Will be `None` if the
            file is not found.
        """

        test = os.path.join(self.confpath, path, conffile)
        files = glob.glob(test, **kwargs)

        n = len(files)
        if n == 0:
            LOGGER.warning("no config file found")
        elif n == 1:
            return files[0]
        else:
            LOGGER.warning('Multiple configuration files match criterion')

    @property
    def h5pyargs(self):
        """
        Property of the HDF5 configurations

        Returns
        -------
        dic : dict
           A dictionary of hte HDF5 config

        """

        return {'compression': self.compression,
                'compression_opts': self.compression_opts}

    def __str__(self):
        """
        Method to override printing
        """
        dic = self.__dict__

        width = 0
        for d in dic.keys():
            width = max(width, len(d))

        s = 'Base Parameters:'
        for k, v in dic.items():
            s += f'\n {k:>{width}} {v}'

        s += '\nConfig Parameters:'

        for k, v in self.items():
            s += f'\n {v}'

        return s


# if __name__=='__main__':
#    x=Config()
#
#    print(x.conffile)
#    print(x.confpath)
#    print(x.fluxscale)
#    x.fluxscale=2
#    print(x.fluxscale)
