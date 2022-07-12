from datetime import datetime


from ...logger import LOGGER
from ...info import __version__,__code__,__email__,__author__

#,__ref__,__refurl__


def add_software_log(hdr):

    hdr.set('CODE',value=__code__,comment='name of software')
    hdr.set('VERSION',value=__version__,comment=f'{__code__} version number')
    hdr.set('AUTHOR',value=__author__,comment=f'author of {__code__}')
            #hdr.set('EMAIL',value=__email__,comment=f'email of author')
    #hdr.set('REF',value=__ref__,
    #        comment=f'publication reference')
    #hdr.set('REFURL',value=__refurl__)
    LOGGER.debug('set the conf path to the output')
    #hdr.set('CONFPATH',value=os.environ[CONFIG_ENV],
    #        comment=f'path to config')
    add_stanza(hdr,"Software Log",before='CODE')


def add_preamble(hdr,filetype):
    now=datetime.now()
    hdr.set("DATE",value=now.strftime("%Y-%m-%d"),
            comment='date this file was written (yyyy-mm-dd)')
    hdr.set("FILETYPE",value=filetype,
            comment='contents of this file')
    add_stanza(hdr,'File Properties',before='DATE')


def add_stanza(hdr,label,**kwargs):
    hdr.set('',value='',**kwargs)
    hdr.set('',value=f'      / {label}',**kwargs)
    hdr.set('',value='',**kwargs)
