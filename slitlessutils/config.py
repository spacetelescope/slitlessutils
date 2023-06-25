
import configparser
import glob
import json
import os
from pathlib import Path
import requests
import socket
import tarfile
from packaging import version
from .core.utilities import headers
from .logger import LOGGER


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

    # the default config filename:
    DEFAULTS = 'defaults.cfg'

    # version filename
    VERSIONFILE = 'version.json'

    # the reference files on box
    REFURL = 'https://stsci.box.com/shared/static/'
    REFFILE = 'fzlb7y36ofi18ziy6mkbyg710stmygjf.gz'

    # local path to store reference files
    REFROOT = os.path.join(Path.home(), '.slitlessutils')

    # enable the singleton
    _instance = None

    def __new__(cls, conffile=None, **kwargs):
        if not cls._instance:

            # make a new config object
            cls._instance = super().__new__(cls)
            cls._instance._refpath = os.getcwd()

            # set the latest reference files
            cls._instance.set_reffiles()

            # load the defaults
            deffile = os.path.join(cls._instance.refpath,
                                   cls._instance.DEFAULTS)

            # get the parser
            cfg = configparser.ConfigParser()

            # read the files
            cfg.read(deffile)

            # if user sets a conffile, then let's use that
            if conffile is not None:
                cfg.read(conffile)

            # set the results to the local dict
            for k in cfg.sections():
                cls._instance[k] = Parameter.from_config(cfg[k])

            # update any keywords
            for k, v in kwargs.items():
                if k in cls._instance:
                    cls._instance[k].value = v

        return cls._instance

    def make_refroot(self):
        """
        Method to make root config directory.
        """

        if not os.path.exists(self.REFROOT):
            os.mkdir(self.REFROOT)

    @property
    def refpath(self):
        return self._refpath

    def set_refpath(self, path):
        if os.path.exists(path):
            LOGGER.info(f'Using reference path: {path}')
            self._refpath = path

            # vers = self.read_versfile()

        else:
            LOGGER.warning(f'Reference path ({path}) does not exist.')

    def read_versfile(self, index=0, filename=None):
        """
        Method to read the version file

        Parameters
        ----------
        index : int, optional
            The index of the version information to read (most recent is 0).
            If `None`, then will return all of the version history.
            Default is 0

        filename : str, optional
            The name of the version file.  If set to `None`, then will use
            the class variable `self.VERSIONFILE`.  Default is None.

        Returns
        -------
        data : dict
            The most recent entry in the version file
        """

        # get the filename
        if filename is None:
            filename = self.VERSIONFILE
        versfile = os.path.join(self.refpath, filename)

        # check if the file exists
        if os.path.exists(versfile):
            # open file as a json
            with open(versfile, 'r') as fp:
                data = json.load(fp)

            # return the requested index
            if index is not None:
                data = data[index]
        else:
            # return empty dict for a non-existent VERSION info
            LOGGER.warning(
                f'Unable to find version file ({versfile}), suspicious behavior may follow.')
            data = {}
        return data

    def set_reffiles(self):
        """
        Set the reference file path in the root configurations
        """

        # grab all the paths from the root directory
        glob_token = os.path.join(self.REFROOT, '*')
        paths = glob.glob(glob_token)

        # check that there are valid paths to parse
        if len(paths) == 0:
            # there are no paths, which means there are no ref files.
            # so we should go fetch some
            if self.retrieve_reffiles():
                # check the paths again
                paths = glob.glob(glob_token)
                if len(paths) == 0:
                    LOGGER.warning("Failed auto setup of reference files")

        # initalize with a nonsense value
        bestvers = version.parse('0.0.0')
        bestpath = None

        # check all the paths and take the highest version
        for path in paths:
            if os.path.isdir(path):
                basepath = os.path.basename(path)
                vers = version.parse(basepath)
                if vers >= bestvers:
                    bestvers = vers
                    bestpath = os.path.join(self.REFROOT, path)
                    
        # set the highest version
        if bestpath is None:
            LOGGER.warning(f'Unable to determine the reference directory in {self.REFROOT}')
        else:
            self._refpath = bestpath

    def retrieve_reffiles(self, refurl=None, reffile=None, update=True):
        """
        Method to retrieve files from box

        Parameters
        ----------
        refurl : str or None, optional
            The URL for the box file.  If value is None, then the
            class variable `Config.REFURL` is used.   Default is None

        reffile : str or None, optional
            The name of the box file.  If value is None, then the
            class variable `Config.REFFILE` is used.  Default is None.

        update : bool, optional
            Flag to update the `refpath` after retrieving.  Default is True

        Returns
        -------
        flag : bool
            Flag if retrieval worked properly

        """

        # first check if internet is alive
        #IP = socket.gethostbyname(socket.gethostname())
        #if IP == '127.0.0.1':
        #    LOGGER.warning(f'Invalid IP: {IP} --- check internect connection')
        #    return False

        # parse the inputs
        if refurl is None:
            refurl = self.REFURL

        if reffile is None:
            reffile = self.REFFILE

        # get the name of the files
        remotefile = refurl+reffile
        localfile = os.path.join(self.REFROOT, reffile)

        # make sure refdir exists
        self.make_refroot()

        # write the remote file into local file
        LOGGER.info(f'Retrieving remote file {remotefile} to {self.REFROOT}')
        req = requests.get(remotefile, timeout=10)
        with open(localfile, 'wb') as f:
            f.write(req.content)


        print(localfile)
        # unpack the local file
        with tarfile.open(localfile) as tarf:
            # security flaw?
            # f.extractall(self.REFROOT)
            print(tarf)
            for f in tarf:
                name = f.name
                relname = os.path.join(self.REFROOT, name)
                absname = os.path.abspath(relname)
                print(absname)
                if absname.startswith(self.REFROOT):
                    print(name)
                    tarf.extract(name, path=self.REFROOT)                    
            
        # the current path
        curpath = os.path.join(self.REFROOT, 'slitlessutils_config')

        # read the version number for this file
        data = self.read_versfile()
        vers = data.get('version','0.0.0')
        
        # rename the directory to the version string
        newpath = os.path.join(self.REFROOT, vers)
        os.rename(curpath, newpath)
        
        # use the new directory?
        if update:
            self.set_refpath(newpath)

        # clean up the downloaded file
        os.remove(localfile)

        return True

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
            # get the value of a configuration parameter
            return self[k].value
        else:
            # return the actual value
            return super(Config, self).__getattribute__(k)

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
        if k == 'refpath':
            # treat the reference file path specially
            self.set_refpath(v)
        elif k in self:
            # set configuration parameters that are dict keywords
            self[k].value = v
        else:
            # deal with normal attribute setting
            super(Config, self).__setattr__(k, v)

    def update_header(self, hdr):
        """
        Method to update a fits header

        Parameters
        ----------
        hdr : `astropy.io.fits.Header`
            The fits header to update
        """

        # put the reference file info into the header
        vers = self.read_versfile()
        hdr.set('refpath', value=str(self.refpath),
                comment='reference file path')
        hdr.set('refvers', value=vers.get('version', 'None'),
                comment='reference file version')
        hdr.set('refdate', value=vers.get('date', 'None'),
                comment='reference file date')
        headers.add_stanza(hdr, 'Reference File Settings', before='refpath')

        # put the configuration info into the header
        first = None
        for k, v in self.items():
            v.update_header(hdr)
            if not first:
                first = v.keyword
        headers.add_stanza(hdr, 'Configuration Settings', before=first)

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

    def get_reffile(self, reffile, path='', **kwargs):
        """
        Method to get a reference file from the local directory tree

        Parameters
        ----------
        reffile : str
            The filename of the reference file

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

        glob_token = os.path.join(self.refpath, path, reffile)
        files = glob.glob(glob_token, **kwargs)

        n = len(files)
        if n == 0:
            LOGGER.warning(f"The requested reference file ({reffile}) is not found.")
        elif n == 1:
            return files[0]
        else:
            LOGGER.warning(f'Multiple reference files found for {reffile}.')

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
                'compression_opts': int(self.compression_opts)}

    def __str__(self):
        """
        Method to override printing
        """

        dat = self.read_versfile()

        # do some stuff for reffiles files
        s = 'Slitlessutils Configuration\n'

        s += '  Reference Files:\n'
        s += f'     refpath: {self.refpath}\n'
        s += f'     version: {dat.get("version")}\n'
        s += f'    use date: {dat.get("date")}\n'

        # get the width of the config parameters
        width = 0
        for d in self.keys():
            width = max(width, len(d))

        # update for the config parameters
        s += '\n  Configuration:'
        for k, v in self.items():
            s += f'\n {k:>{width}}: {v.value}'

        # add new line, just because I like it
        s += '\n'

        return s


# if __name__=='__main__':
#    x=Config()
#
#    print(x.conffile)
#    print(x.confpath)
#    print(x.fluxscale)
#    x.fluxscale=2
#    print(x.fluxscale)
