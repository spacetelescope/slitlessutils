from pathlib import Path

from astropy.io import fits
from drizzlepac import astrodrizzle



def _get_instrument_defaults(file, instrument=None):
    # Args taken from INS HSTaXe FullFrame Cookbooks
    instrument_args = {
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

    # Immediately return if the instrument was provided
    if instrument in instrument_args.keys():
        return instrument_args[instrument]

    with fits.open(file, mode='readonly') as hdul:
        h = hdul[0].header
        telescope = h['TELESCOP']
        if telescope == 'HST':
            instrument = h['INSTRUME']
            if instrument == 'WFC3':
                csmid = h['CSMID']
                if csmid == 'IR':
                    return instrument_args['ir']
                elif csmid == 'UVIS':
                    return instrument_args['uvis']
                else: 
                    raise ValueError(f'Invalid CSMID: {csmid}')
            elif instrument == 'ACS':
                if h['DETECTOR'] == 'WFC':
                    return instrument_args['acs']
                elif h['DETECTOR'] == 'SBC':
                    raise ValueError(f'SBC does not need drizzling')
            else:
                raise ValueError(f'\'{instrument}\' is not supported')
        else: 
            raise ValueError(f'\'{telescope}\' is not supported')


def drizzle(files, instrument=None, outdir=Path().absolute(), **kwargs):
    '''
    Runs AstroDrizzle as part of Cosmic Ray handling. Passes any additional arguments
    not listed here directly to AstroDrizzle. Any user-specified arguments will override
    

    Parameters
    ----------
    files : str or list of str
        The list of files to be processed by AstroDrizzle
    
    instrument : {'ir', 'uvis', 'acs'}, optional
        One of three instruments on HST of which to apply a set of default arguments
        captured from the official HSTaXe FullFrame Cookbooks. If not provided, will attempt
        to be automatically detected from the header of the FIRST input file

    outdir : str or `pathlib.Path`, optional
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
    if isinstance(files, list):
        file_to_check = files[0]
    elif isinstance(files, str):
        with open(files, 'r') as f:
            file_to_check = f.readline()
    drizzle_kwargs.update(_get_instrument_defaults(file_to_check, instrument))
    # Finally override any args with the ones the user supplied
    drizzle_kwargs.update(kwargs)
    # Prepend outdir to output
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    drizzle_kwargs['output'] = str(outdir / drizzle_kwargs['output'])

    return astrodrizzle.AstroDrizzle(files, **drizzle_kwargs)
