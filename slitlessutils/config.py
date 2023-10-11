import glob
import json
import os
import shutil
from pathlib import Path
from astropy.utils.data import download_file
import tarfile
from packaging import version
import textwrap


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
    def from_dict(cls, p):
        """
        Classmethod to load config parameter from a config file

        Parameters
        ----------
        p : dict (or dict like)
            A dictionary to load the parameter from

        Returns
        -------
        obj : `Parameter`
            The output parameter
        """

        return cls(p['keyword'], p['value'], p['comment'])

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
    DEFAULTS = 'defaults.json'

    # version filename
    VERSIONFILE = 'version.json'

    # the reference files on box
    # REFURL = 'https://stsci.box.com/shared/static/'
    # REFFILE = 'fzlb7y36ofi18ziy6mkbyg710stmygjf.gz'
    REFURL = 'https://data.science.stsci.edu/redirect/slitlessutils/'
    REFDB = 'slitlessutils_config_versions.json'

    # local path to store reference files
    REFROOT = os.path.join(Path.home(), '.slitlessutils')

    # flag to timeout file downloading in seconds
    TIMEOUT = 15

    # enable the singleton
    _instance = None

    # a basic init method.  Seems necessary to implement kwargs
    def __init__(self, *args, **kwargs):
        pass

    # the main
    def __new__(cls, conffile=None, **kwargs):

        if not cls._instance:

            # make a new config object
            cls._instance = super(Config, cls).__new__(cls)
            cls._instance._refpath = os.getcwd()

            # first check if the refroot exists
            cls._instance.make_refroot()

            # set the latest reference files
            _ = cls._instance.set_reffiles()

            # load the defaults
            deffile = os.path.join(cls._instance.refpath,
                                   cls._instance.DEFAULTS)

            # read the defaults as json
            if os.path.exists(deffile):
                with open(deffile, 'r') as f:
                    cfg = json.load(f)
            else:
                cfg = {}

            # if a conffile is given, read and update with that
            if conffile:
                with open(conffile, 'r') as f:
                    custom = json.load(f)
                cfg.update(custom)

            # load each keyword into the self
            cls._instance._keys = tuple(cfg.keys())
            for k, v in cfg.items():
                cls._instance[k] = Parameter.from_dict(v)

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
    def refversion(self):
        data = self.read_versfile()
        return data.get('version')

    @property
    def refdate(self):
        data = self.read_versfile()
        return data.get('date')

    @property
    def refpath(self):
        return self._refpath

    def set_refpath(self, path):
        """
        Method to set the path to the reference files

        Parameters
        ----------
        path : str
            The path to the reference files.  Code will check that the
            path exists.  It is not required, but code will look for a
            `version.json` file in that directory

        """

        if not path.startswith(self.REFROOT):
            path = os.path.join(self.REFROOT, path)

        if os.path.exists(path):
            # ensure that the last char is a separator
            if path[-1] != os.sep:
                path += os.sep

            # report result and set the variable
            LOGGER.info(f'Using reference path: {path}')
            self._refpath = path

        else:
            LOGGER.warning(f'Reference path ({path}) does not exist.')

    def read_versfile(self, filename=None):
        """
        Method to read the version file

        Parameters
        ----------
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
            filename = os.path.join(self.refpath, self.VERSIONFILE)

        # check if the file exists
        if os.path.exists(filename):
            # open file as a json
            with open(filename, 'r') as fp:
                data = json.load(fp)

        else:
            # return empty dict for a non-existent VERSION info
            msg = f'Unable to find version file ({filename}).  This should be concerning.'
            LOGGER.warning(msg)
            data = {}

        return data

    def set_reffiles(self, refversion=None):
        """
        Set the reference file path in the root configurations

        Parameters
        ----------
        refversion : str, `version.Version`, optional
            Name of the version to use.  If `None`, then use latest.
            default is None.
        """

        # for messaging:
        failure = "Cannot set refpath:"

        # process based on the data type of the refversion
        if refversion is None:
            # if refvserion is `None`, then compute the latest-greatest.

            # initalize with a nonsense value
            bestvers = version.parse('0.0.0')
            bestpath = None

            # get the paths in the root dir
            glob_token = os.path.join(self.REFROOT, '*'+os.sep)
            paths = glob.glob(glob_token)

            # check for empty paths
            if len(paths) == 0:
                msg = f"{failure} there are no reference files in " + \
                    "{self.REFROOT}, trying to fetch some."
                LOGGER.warning(msg)

                # try retrieving some files
                reffile = self.retrieve_reffiles(update=True)
                if not reffile:
                    LOGGER.error("Cannot set any reffiles")
            else:
                # check all the paths and take the highest version
                for path in paths:
                    basepath = os.path.basename(path[:-1])
                    try:
                        vers = version.parse(basepath)
                    except version.InvalidVersion:
                        vers = None
                    if vers and (vers >= bestvers):
                        bestvers = vers
                        bestpath = os.path.join(self.REFROOT, path)

                # ok now set the best path
                if bestpath is None:
                    LOGGER.warning(f"{failure} latest version not found in {self.REFROOT}")
                    return
                else:
                    self.set_refpath(bestpath)

        else:
            # if not `None`, then sort it out
            if isinstance(refversion, str):
                # if a string, then use it as is
                vers = refversion
            elif isinstance(refversion, version.Version):
                # if a Python version object, then distill as a string
                vers = str(refversion)
            else:
                # invalid data type
                LOGGER.warning(f"{failure} invalid datatype for refversion ({refversion})")
                return

            # now look for this refversion in the archive of reference files:
            path = os.path.join(self.REFROOT, refversion)
            if os.path.isdir(path):
                # valid path was found.  Use it
                self.set_refpath(path)
            else:
                # no valid path was found.  give an error
                LOGGER.warning(f"{failure} refversion not found ({refversion})")

        # return a value to indicate success
        return self._refpath

    def help_refmanifest(self):
        """
        Method to print properties of the reference-file manifest to
        the screen

        """
        # indentations
        used = ' * '
        others = '   '

        # get manifest and print item by item
        man = self.retrieve_refmanifest()
        if not man:
            return

        ref_version_found = False
        print('Reference File Manifest:\n')
        for vers, data in man.items():

            # get the version number, indentation, and format them
            vers = str(vers)
            if vers == self.refversion:
                indent = used
                ref_version_found = True
            else:
                indent = others

            # print the main line
            print(f'{indent}{vers:<9} {data["file"]:<40} {data["date"]}')

            # print the notes field
            for text in textwrap.wrap(data["notes"]):
                print(f'{others}{text}')
            print()

        if not ref_version_found:
            LOGGER.warning(
                f"Manifest does not include current reference version: {self.refversion}")

    def retrieve_refmanifest(self, return_latest=False):
        """
        Method to read the manifest file about all the reference files

        Parameters
        ----------
        return_latest : bool, optional
            Flag to also return the highest (ie. latest) version number.
            Default is False

        Returns
        -------
        manifest : dict
            A dictionary of the possible reference files.  Keywords are
            version numbers, values are the filenames.

        version : `version.Version`, optional
            The version number of the highest (ie. latest) version, if
            requested.

        """

        # download the manifest file
        try:
            f = download_file(self.REFURL+self.REFDB, timeout=self.TIMEOUT,
                              show_progress=False)
        except TimeoutError:
            LOGGER.warning(f'Retrieving manifest timed out in {self.TIMEOUT} s.')
            return

        # open the manifest as a json file
        with open(f, 'r') as fp:
            data = json.load(fp)

        # output variables:
        manifest = {}
        latest = version.parse('0.0.0')

        # parse the manifest
        for datum in data:
            vers = version.parse(datum.get('version', '0.0.0'))
            if vers >= latest:
                latest = vers

            # remove the version variable (since it's the keyword) and
            # save to manifest dict for output
            if 'version' in datum:
                del datum['version']
            manifest[vers] = datum

        # sort out what to return
        if return_latest:
            return manifest, latest
        else:
            return manifest

    def retrieve_reffiles(self, refversion=None, update=True):
        """
        Method to retrieve files from box

        Parameters
        ----------
        refversion : str or `version.Version`, optional
            The name of the reference file version info.  If `None`, then the
            latest (ie. highest) version are used.  Default is None.

        update : bool, optional
            Flag to update the `refpath` after retrieving.  Default is True

        Returns
        -------
        reffile : str
            name of the reference library downloaded or `None` if retrieval failed.

        """

        # read the manifest
        man, latest = self.retrieve_refmanifest(return_latest=True)
        if refversion is None:
            reffile = man[latest]['file']
        else:
            if isinstance(refversion, str):
                refversion = version.parse(refversion)
            elif isinstance(refversion, version.Version):
                pass
            else:
                LOGGER.warning(f"Unsupported datatype for refversion ({refversion})")
                return

            if refversion in man:
                reffile = man[refversion]['file']
            else:
                LOGGER.warning(f"Refversion ({refversion}) not found in manifest.")
                return

        # get the name of the files
        remotefile = self.REFURL+reffile
        localfile = os.path.join(self.REFROOT, reffile)

        # write the remote file into local file
        LOGGER.info(f'Retrieving remote file {remotefile} to {self.REFROOT}')
        try:
            tmpfile = download_file(remotefile, timeout=self.TIMEOUT)
        except BaseException:
            LOGGER.warning(f"Unable to retrieve server file ({remotefile})")
            return
        except TimeoutError:
            LOGGER.warning(f'Download timed-out in {self.TIMEOUT} s.')
            return

        # move the file from a tmp dir
        shutil.move(tmpfile, localfile)

        # unpack the local file
        with tarfile.open(localfile) as tarf:
            for f in tarf:
                name = f.name
                relname = os.path.join(self.REFROOT, name)
                absname = os.path.abspath(relname)
                if absname.startswith(self.REFROOT):
                    tarf.extract(name, path=self.REFROOT)

        # delete the local file now that we're done with it
        os.remove(localfile)

        # the current path
        curpath = os.path.join(self.REFROOT, 'slitlessutils_config')

        # read the version number for this file
        versfile = os.path.join(curpath, self.VERSIONFILE)
        versdata = self.read_versfile(filename=versfile)
        vers = versdata.get('version', '0.0.0')

        # get name of the new directory
        newpath = os.path.join(self.REFROOT, vers)

        # if new path is full, then empty it
        if os.path.isdir(newpath) and os.listdir(newpath):
            shutil.rmtree(newpath)

        # rename path
        os.rename(curpath, newpath)

        # use the new directory?
        if update:
            self.set_refpath(newpath)

        return remotefile

    def __setitem__(self, k, v):
        """
        Method to override the dict-like access

        Parameters
        ----------
        k : str
            Name of the dict keyword

        v : any type
            The value of the attribute
        """

        if isinstance(v, Parameter):
            super().__setitem__(k, v)
        elif k in self._keys:
            self[k].value = v
        else:
            pass

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
            try:
                # return the actual value
                return super(Config, self).__getattribute__(k)
            except AttributeError:
                LOGGER.warning(f"Attribute {k} not found.")
                return

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
        Method to write config file to disk as a *.json type

        Parameters
        ----------
        filename : str
            The name of the config file

        """

        ext = os.path.splitext(filename)[-1][1:]
        if ext in ('json',):
            data = {k: dict(v) for k, v in self.items()}

            with open(filename, 'w') as f:
                json.dump(data, f, indent=4)

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
                'compression_opts': self.compression_opts}

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
