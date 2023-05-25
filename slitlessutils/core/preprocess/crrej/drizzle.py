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
