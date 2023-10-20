import logging
import os
import sys
from datetime import datetime
from multiprocessing import current_process

from .info import __code__

LOGGERNAME = __code__
KNOWNISSUE = 15
RTD = 'https://slitlessutils.readthedocs.io/en/latest/knownissues.html'


class BaseFormatter(logging.Formatter):
    """
    Base class to implement custom formatting
    """

    DEFAULT = '%(message)s'

    def format(self, record):
        """
        Mandatory method to do the formatting

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

class STDOUTFormatter(BaseFormatter):
    """
    Class to log to stdout, which enables fancy color printing
    """

    FORMATS = {logging.NOTSET: '%(message)s',
               logging.DEBUG: ("\033[34;1;3mDEBUG (%(module)s:%(lineno)s)> "
                               "\033[00m\033[34;3m%(message)s\033[00m"),
                logging.INFO: "\033[32;1mINFO %(processName)s> \033[00m\033[32m%(message)s\033[00m",
                KNOWNISSUE: f"\n\033[36;1mKNOWN ISSUE (%(module)s:%(lineno)s)> \033[00m\033[36m%(message)s\nSee: \033[04m{RTD}\033[00m\n",
                logging.WARNING: ("\033[93;1mWARNING %(processName)s> "
                                 "\033[00m\033[93m%(message)s\033[00m"),
                logging.ERROR: ("\033[91;5;1mERROR %(processName)s>\033[00m"
                              "\033[91m %(message)s\033[00m"),
                logging.CRITICAL: "\033[91;1;5;7mCRITICAL %(processName)s> %(message)s\033[00m"}

class FileFormatter(BaseFormatter):
    """
    Class to enable ascii-file printing, which has extra information.
    """

    FORMATS = {logging.NOTSET: "%(message)s",
                logging.DEBUG: "   DEBUG (%(module)s:%(lineno)s)> %(message)s",
                logging.INFO: "    INFO> %(message)s",
                KNOWNISSUE: f"\nKNOWNISSUE (%(module)s:%(lineno)s)> %(message)s\nSee: {RTD}\n",
                logging.WARNING: " WARNING> %(message)s",
                logging.ERROR: "   ERROR> %(message)s",
                logging.CRITICAL: "CRITICAL> %(message)s"}


def setLevel(level, logger=None):
    """
    Helper function to set the logging level of all the handlers
    
    Parameters
    ----------
    level : int
        The level to set

    logger : `logging.Logger()`, optional
        A logger to set the level for, if None, then use the __code__ variable.
        Default is None.    
    """

    if logger is None:
        logger = logging.getLogger(LOGGERNAME)
    
    logger.setLevel(level)
    for handler in logger.handlers:
        handler.setLevel(level)


def initialize(level=10, asciilog=True, stdout=True, stderr=False):
    """
    Initialize the custom Logger

    Parameters
    ----------
    level : int or logging level name, optional.
        The default level to log.  Default = 10

    asciilog : bool, optional
        Flag to log to an ascii file named "slitlessutils.log".  Default is True

    stdout : bool, optional
        Flag to log to STDOut.  Default is True

    stderr : bool, optional
        Flag to log to STDErr.  Default is False

    Returns
    -------
    log : `logging.Logger()`
        The initialized logger.

    """

    # grab the logger, and check if it already exists
    log = logging.getLogger(LOGGERNAME)
    if log.hasHandlers():
        return log

    # make a function for known issues
    log.knownissue = lambda msg, *args, **kwargs: log._log(KNOWNISSUE, msg, args, **kwargs) if log.isEnabledFor(KNOWNISSUE) else None
                
    # create a file handler
    if asciilog:

        # sort out the log filename
        procname = current_process().name
        if procname == 'MainProcess':
            logfile = f'{LOGGERNAME}.log'
        else:
            logfile = f'{LOGGERNAME}-{procname}.log'
            
        # zap a current file
        if os.path.exists(logfile):
            os.remove(logfile)

        # establish the handler
        filehandler = logging.FileHandler(logfile)
        filehandler.setFormatter(FileFormatter())
        log.addHandler(filehandler)

    # create a stdout handler
    if stdout:
        stdouthandler = logging.StreamHandler(sys.stdout)
        stdouthandler.setFormatter(STDOUTFormatter())
        log.addHandler(stdouthandler)
        
    # create a stderr handler
    if stderr:
        stderrhandler = logging.StreamHandler(sys.stderr)
        stderrhandler.setFormatter(STDOUTFormatter())
        log.addHandler(stderrhandler)

    # set the levels
    setLevel(level, logger=log)

    # do some pre-logging
    log.info(f"Logger ({LOGGERNAME}) with {level=} at {datetime.now()}")

    return log

LOGGER = initialize()

#if __name__=='__main__':        
#    
#    LOGGER = initialize()
#    q=LOGGER.info('test')
#    a=LOGGER.knownissue('known issue')
#    
#    print(q,a)
