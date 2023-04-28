import logging
import sys
import os
from datetime import datetime
from multiprocessing import current_process

from .info import __code__

"""
A collection of functions and classes to implement the logging.

The typical use case is
>>> import logger
>>> logger.start_logging()
>>> logger.LOGGER.warning('a message')

or after the logger has been initialized
>>> from logger import LOGGER
>>> LOGGER.info('this is an informational message')

But one can turn off the logging:
>>> import logger
>>> logger.end_logging()

Or changing the logging level
>>> import logger
>>> logger.setLevel(20)

"""


# the default logger
LOGGER = logging.getLogger(__code__)


class MyFormatter(logging.Formatter):
    """
    Base class to implement custom formatting
    """

    DEFAULT = '%(message)s'

    def format(self, record):
        """
        Mandatory method to do the formating

        Parameters
        ----------
        record : `logging.LogRecord`
            the content to format

        Returns
        -------
        fmt : `logging.LogRecord`
            the formatted record
        """

        self._style._fmt = self.FORMATS.get(record.levelno, self.DEFAULT)
        fmt = logging.Formatter.format(self, record)
        return fmt


# %(processName)s
# %(process)d

class STDOUTFormatter(MyFormatter):
    """
    Class to log to stdout, which enables fancy color printing
    """

    FORMATS = {logging.NOTSET: '%(message)s',
               logging.DEBUG: "\033[34;1;3mDEBUG (%(module)s:%(lineno)s)> \033[00m\033[34;3m%(message)s\033[00m",
               logging.INFO: "\033[32;1mINFO %(processName)s> \033[00m\033[32m%(message)s\033[00m",
               logging.WARNING: "\033[93;1mWARNING %(processName)s> \033[00m\033[93m%(message)s\033[00m",
               logging.ERROR: "\033[91;5;1mERROR%(processName)s>\033[00m\033[91m %(message)s\033[00m",
               # logging.ERROR :"\033[91;1mERROR %(processName)s>
               # \033[00m\033[91;5m%(message)s\033[00m",
               logging.CRITICAL: "\033[91;1;5;7mCRITICAL %(processName)s> %(message)s\033[00m"}


class FileFormatter(MyFormatter):
    """
    Class to enable ascii-file printing, which has extra information.
    """

    FORMATS = {logging.NOTSET: "%(message)s",
               logging.DEBUG: "   DEBUG (%(module)s:%(lineno)s)> %(message)s",
               logging.INFO: "    INFO> %(message)s",
               logging.WARNING: " WARNING> %(message)s",
               logging.ERROR: "   ERROR> %(message)s",
               logging.CRITICAL: "CRITICAL> %(message)s"}


def setLevel(level):
    """
    Helper function to set the logging level of all the handlers

    Parameters
    ----------
    level : int
       The level to set
    """

    log = logging.getLogger(__code__)
    log.setLevel(level)
    for handler in log.handlers:
        handler.setLevel(level)


def start_logging(asciilog=True, stdout=True, stderr=False, level=10):
    """
    Helper function to start the logging

    Parameters
    ----------
    asciilog : bool, optional
       Flag to log to an ascii file.  Default is True

    stdout : bool, optional
       Flag to log to the stdout.  Default is True

    stderr : bool, optional
       Flag to log to the stderr.  Default is False

    level : int, optional
       Level to initialize for logging.  Default is 10

    Notes
    -----
    If logging to an ascii file, then the filename will be automatically
    generated.  If the process begins in the 'main', then it will be simply:
    `slitlessutils.log`.  If the process begins in a child process (say
    from `multiprocessing`), then the file will be `slitlessutils-NUMBER.log`
    Obviously, these files will be overwritten each time, as they are clearly
    not unique.

    """

    # dct=logging.Logger.manager.loggerDict
    # if __code__ in dct and :
    #    print(dct[__code__])
    # else:
    #    print(dct.values())
    #
    if LOGGER.level == level:
        return

    # sort out the log file
    procname = current_process().name
    if procname == 'MainProcess':
        logfile = f'{__code__}.log'
    else:
        logfile = f'{__code__}-{procname}.log'
    if os.path.exists(logfile):
        os.remove(logfile)

    # get the logger
    log = logging.getLogger(__code__)

    # create a file handler
    if asciilog:
        filehandler = logging.FileHandler(logfile)
        filehandler.setLevel(level)
        filehandler.setFormatter(FileFormatter())
        log.addHandler(filehandler)

    # create a stdout handler
    if stdout:
        stdouthandler = logging.StreamHandler(sys.stdout)
        stdouthandler.setLevel(level)
        stdouthandler.setFormatter(STDOUTFormatter())
        log.addHandler(stdouthandler)

    # create a stderr handler
    if stderr:
        stderrhandler = logging.StreamHandler(sys.stderr)
        stderrhandler.setLevel(level)
        stderrhandler.setFormatter(STDOUTFormatter())
        log.addHandler(stderrhandler)

    # change teh default logging level
    setLevel(level)

    # do some pre-logging
    log.info(f"Starting logger with {level=}")
    log.info(f'Importing {__code__} at {datetime.now()}')


def end_logging():
    """
    Helper function to discontinue the logging
    """

    # get the logger
    log = logging.getLogger(__code__)
    log.handlers.clear()
    del log


def initialize_logger(level=10, start=False, **kwargs):
    """
    Helper function to initialize the loggers

    Parameters
    ----------
    level : int, optional
        The initial level of the logger.  Default is 10.

    start : bool, optional
        Flag to additionally start the logging instead of only initializing
        the logging handlers.  Default is False

    kwargs : dict, optional
        Dictionary to pass to `start_logging()`

    """

    # create a logger
    log = logging.getLogger(__code__)

    if start:
        start_logging(level=level, **kwargs)
    else:
        # null log handler
        nullhandler = logging.NullHandler()
        nullhandler.setLevel(level)
        log.addHandler(nullhandler)

    return log


if __name__ == '__main__':

    logger = initialize_logger('slitless', level=30)

    logger.info('info message')
    setLevel('slitless', 1)

    logger.warning('info message\n second line \n third line')

    # to capture print message via logging
    # from contextlib import redirect_stdout
    # logger.write = lambda msg: logger.log(11,msg) if msg != '\n' else None
    #
    # with redirect_stdout(logger):
    #    print('this is a test')
