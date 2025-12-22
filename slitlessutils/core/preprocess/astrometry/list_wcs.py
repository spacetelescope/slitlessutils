import argparse as ap
from glob import glob

from astropy.io import fits


def list_wcs(arg):
    '''
    Method to list the WCSs present in a host of flt/flc files by
    printing to the screen.  The live WCSNAME will be shown in boldface.

    Parameters
    ----------
    arg : str
       Token to pass to `glob.glob()` to find all files.

    Returns
    -------
    None.
    '''

    files = glob(arg)

    for f in files:

        with fits.open(f, mode='readonly') as hdul:
            h0 = hdul[0].header

            obs = (h0['TELESCOP'], h0['INSTRUME'], h0['DETECTOR'])
            if obs == ('HST', 'ACS', 'WFC'):
                filtername = h0['FILTER1']
            elif obs == ('HST', 'ACS', 'SBC'):
                filtername = h0['FILTER1']
            elif obs == ('HST', 'WFC3', 'UVIS'):
                filtername = h0['FILTER']
            elif obs == ('HST', 'WFC3', 'IR'):
                filtername = h0['FILTER']
            else:
                raise ValueError(f"Cannot find observation {obs}")

            # find all possible WCS
            data = {}
            names = []
            for hdu in hdul:
                hdr = hdu.header
                extname = hdr.get("EXTNAME")
                extver = hdr.get("EXTVER")

                if extname == 'SCI' and isinstance(extver, int):
                    data[extver] = {}
                    for k, v in hdr.items():
                        if k.startswith("WCSNAME"):
                            data[extver][k] = v
                            names.append(k)
            names = set(names)

            # print a tabulated view
            print(f' {f:43} {filtername}')
            print('+' + '-' * 10 + '+' + '-' * 32 + '+' + '-' * 32 + '+')

            vers = ' | '.join(f'{k:30}' for k in data.keys())

            print(f'| WCSVERS  | {vers} |')
            print('+' + '-' * 10 + '+' + '-' * 32 + '+' + '-' * 32 + '+')
            for name in sorted(names):

                bold = (name == 'WCSNAME')
                if bold:
                    row = f'|\033[1m {name:8} \033[0m| '
                else:
                    row = f'| {name:8} | '

                for extver in data.keys():
                    wcs = data[extver].get(name, '')
                    if bold:
                        row += f'\033[1m{wcs:>30}\033[0m | '
                    else:
                        row += f'{wcs:>30} | '

                # show the row
                print(row)

            print('+' + '-' * 10 + '+' + '-' * 32 + '+' + '-' * 32 + '+')
            print('\n')


def main():

    # set up the argument
    p = ap.ArgumentParser(description='List the WCSs for a collection of files')

    # add suffix to the list
    p.add_argument("suffix", help='The fits file suffix, example: "flc"')

    # read the arguments
    args = p.parse_args()

    # call the function
    list_wcs(args.suffix)
