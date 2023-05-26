from drizzlepac import astrodrizzle
from pathlib import Path

# Args taken from INS HSTaXe FullFrame Cookbooks
COMMON_ARGS = {
    'ir': {
        'mdriztab': True,
        'preserve': False,
        'skysub': False,
        'final_fillval': None
    },
    'uvis': {
        'mdriztab': True,
        'preserve': False,
        'skysub': True,
        'skymethod': 'match',
        'final_fillval': None
    },
    'acs': {
        'build': True,
        'mdriztab': True,
        'in_memory': False,
        'skysub': False,
        'driz_separate': False,
        'median': False,
        'blot': False,
        'driz_cr': False,
        'driz_sep_wcs': False,
        'driz_combine': True,
        'final_wcs': False
    }
}


def drizzle(files, instrument=None, outdir=Path().absolute(), **kwargs):
    '''
    Runs AstroDrizzle as part of Cosmic Ray handling. Passes any additional arguments
    not listed here directly to AstroDrizzle. Any user-specified arguments will override
    

    Parameters
    ----------
    files : str to a file catalog, or a python list of files
        The list of file to be processed by AstroDrizzle
    
    instrument : str (one of 'ir', 'uvis', 'acs')
        One of three instruments on HST of which to apply a set of default arguments
        captured from the official HSTaXe FullFrame Cookbooks

    outdir : str or `pathlib.Path`
        A directory to write the final rectified mosaics to.
        By default, the current working directory
    '''
    # Start with known defaults
    drizzle_kwargs = {
        'output': 'su_drizzle',
        'overwrite': False,
        'preserve': False,
        'clean': True
    }
    # Apply instrument-specific defaults
    drizzle_kwargs.update(COMMON_ARGS.get(instrument, {}))
    # Finally override any args with the ones the user supplied
    drizzle_kwargs.update(kwargs)
    # Prepend outdir to output
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    drizzle_kwargs['output'] = str(outdir / drizzle_kwargs['output'])

    astrodrizzle.AstroDrizzle(files, **drizzle_kwargs)
