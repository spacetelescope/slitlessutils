import os

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

import slitlessutils as su

ROOT = 'pointsources'

# parameters for the set up of the instrument
TELESCOPE = 'HST'
INSTRUMENT = 'WFC3'
DETECTOR = 'IR'
WRANGE = (7500., 11000.)

SUFFIX = 'flt'
SEDPATH = 'cdbs'


def getPrefix(ins):
    if ins == 'WFC3':
        return 'i'
    elif ins == 'ACS':
        return 'j'
    else:
        raise NotImplementedError("Invalid instrument")


def extract_multi():
    ins = getPrefix(INSTRUMENT)

    # load the grism images
    data = su.wfss.data.WFSSCollection.from_glob(f'{ins}*_{SUFFIX}.fits.gz')

    # load the sources into SU
    sources = su.sources.SourceCollection(f'{ROOT}_seg.fits',
                                          f'{ROOT}_sci.fits',
                                          local_back=False)

    # test grouping
    grp = su.modules.Group(ncpu=1, orders=('+1',))
    groups = grp(data, sources)
    groups.plot(f'{ROOT}_grp.pdf')

    # run the single-file extraction
    ext = su.modules.Multi('+1', (-5., 0, 0.1), root=ROOT, algorithm='grid')
    res = ext(data, sources, groups=groups)
    return res


def makeplot():
    fluxscale = su.config.Config().fluxscale

    AA = '\\mathrm{\\AA}'
    filename = f'{ROOT}_x1d.fits'
    with fits.open(filename, mode='readonly') as hdul:

        for i, hdu in enumerate(hdul[1:]):
            data = hdu.data

            nnan = np.sum(np.isnan(data['flam']))
            ntot = len(data['flam'])

            if nnan != ntot:
                segid = hdu.header['SEGID']

                xerr = (data['lamb'][1] - data['lamb'][0]) / 2
                xerr = np.full_like(data['lamb'], xerr)

                plt.plot(data['lamb'], data['flam'], color='darkgrey')
                plt.errorbar(data['lamb'], data['flam'],
                             xerr=xerr, yerr=data['func'],
                             color='black', label='multi_extract', fmt='o')

                specfile = os.path.join('spectra', f'{segid}_1.csv')
                l, f = np.loadtxt(specfile, usecols=(0, 1), delimiter=',',
                                  skiprows=1, unpack=True)

                g = np.where((WRANGE[0] <= l) & (l <= WRANGE[1]))
                l = l[g]
                f = f[g] / fluxscale
                plt.plot(l, f, color='red', label='CRDS')

                ymin = np.amin(f)
                ymax = np.amax(f)
                plt.ylim((0.95 * ymin, 1.05 * ymax))
                plt.ylabel(r'$f_\lambda$ ($10^{-17}$ erg/s/cm$^2$/' + AA + ')')

                plt.xlim(*WRANGE)
                plt.xlabel(r'$\lambda_\mathrm{obs}$ (' + AA + ')')

                plt.legend()
                plt.tight_layout()
                plt.show()


def run_all(plot=False):

    extract_multi()

    if plot:
        makeplot()


if __name__ == '__main__':
    run_all()
