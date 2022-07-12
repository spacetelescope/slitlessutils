import logging
import sys
import os
from datetime import datetime
from multiprocessing import current_process

from .info import __code__


LOGGER=logging.getLogger(__code__)




class MyFormatter(logging.Formatter):
    DEFAULT='%(message)s'
    
    def format(self, record):
        self._style._fmt = self.FORMATS.get(record.levelno,self.DEFAULT)
        return logging.Formatter.format(self, record)
    

# %(processName)s
# %(process)d

class STDOUTFormatter(MyFormatter):
    FORMATS = {logging.NOTSET: '%(message)s',
               logging.DEBUG :"\033[34;1;3mDEBUG (%(module)s:%(lineno)s)> \033[00m\033[34;3m%(message)s\033[00m",
               logging.INFO :"\033[32;1mINFO %(processName)s> \033[00m\033[32m%(message)s\033[00m",
               logging.WARNING :"\033[93;1mWARNING %(processName)s> \033[00m\033[93m%(message)s\033[00m",
               logging.ERROR :"\033[91;1mERROR %(processName)s> \033[00m\033[91;5m%(message)s\033[00m",
               logging.CRITICAL :"\033[91;1;5;7mCRITICAL %(processName)s> %(message)s\033[00m"}






class FileFormatter(MyFormatter):
    FORMATS = {logging.NOTSET: "%(message)s",
               logging.DEBUG   :"   DEBUG (%(module)s:%(lineno)s)> %(message)s",
               logging.INFO    :"    INFO> %(message)s",
               logging.WARNING :" WARNING> %(message)s",
               logging.ERROR   :"   ERROR> %(message)s",
               logging.CRITICAL:"CRITICAL> %(message)s"}



def initialize_logger(name,level=10):

    
    
    # create a logger
    log=logging.getLogger(name)
    


    # sort out the log file
    procname=current_process().name
    if procname == 'MainProcess':
        logfile=f'{name}.log'
    else:
        logfile=f'{name}-{procname}.log'
    if os.path.exists(logfile):
        os.remove(logfile)
        
    # create a file handler
    filehandler=logging.FileHandler(logfile)
    filehandler.setLevel(level)
    filehandler.setFormatter(FileFormatter())
    log.addHandler(filehandler)

    # create a stdout handler
    stdouthandler=logging.StreamHandler(sys.stdout)
    stdouthandler.setLevel(level)
    stdouthandler.setFormatter(STDOUTFormatter())
    log.addHandler(stdouthandler)

    
    # log to stderr
    #stderrhandler=logging.StreamHandler(sys.stderr)
    #stderrhandler.setLevel(level)
    #stderrhandler.setFormatter(STDOUTFormatter())
    #log.addHandler(stderrhandler)

    # change teh default logging level
    #log.setLevel(level)
    setLevel(name,level)
    
    # do some pre-logging
    log.info(f'Importing {__code__} at {datetime.now()}')

    
    return log


def setLevel(name,level):
    log=logging.getLogger(name)
    log.setLevel(level)
    for handler in log.handlers:
        handler.setLevel(level)


if __name__=='__main__':

    logger=initialize_logger('slitless',level=30)

    logger.info('info message')
    setLevel('slitless',1)
    
    logger.warning('info message\n second line \n third line')




    ## to capture print message via logging
    #from contextlib import redirect_stdout
    #logger.write = lambda msg: logger.log(11,msg) if msg != '\n' else None
    #
    #with redirect_stdout(logger):
    #    print('this is a test')
