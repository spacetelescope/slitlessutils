import slitlessutils as su
from glob import glob

from .parameters import ROOT

def extract_single(**kwargs):
    
    # load data into SU
    data=su.data.wfss.WFSSCollection.from_list(glob(f'*{ROOT}_flc.fits.gz'))

    # load the sources into SU
    sources=su.sources.SourceCollection(f'{ROOT}_seg.fits',
                                        f'{ROOT}_img.fits')

    
    # now run the simulator
    ext = su.modules.Single(**kwargs)
    res = ext(data,sources)

    return res

