import slitlessutils as su

from .parameters import ROOT

def simulate(**kwargs):
    

    # load data into SU
    data=su.data.wfss.WFSSCollection.from_wcsfile(f'{ROOT}_wcs.csv')

    # load the sources into SU
    sources=su.sources.SourceCollection(f'{ROOT}_seg.fits',
                                        f'{ROOT}_img.fits',
                                        sedfile=f'{ROOT}_seds.fits')
    sources.write_seds(path='seds')

    
    # tabulate everything first!
    tab=su.modules.Tabulate(**kwargs)
    pdtfiles=tab(data,sources)

    # now run the simulator
    sim = su.modules.Simulate(**kwargs)
    imgfiles = sim(data,sources)

    return imgfiles



