import os
import warnings

from astropy.io import fits
from astropy.modeling import models
from astropy.wcs import WCS
import numpy as np

import slitlessutils as su

ROOT = 'pointsources'

# parameters for the set up of the instrument
TELESCOPE = 'HST'
INSTRUMENT = 'WFC3'
DETECTOR = 'IR'
DISPERSER = 'G102'
BLOCKING = ''
FILTER = 'F105W'
WRANGE = (7500., 11000.)
SUFFIX = 'flt'
EXPTIME = 1000.

# Notional Field Parameters
RA = 53.
DEC = -27.
PIXSCL = 0.06
SIZE = (100, 100)
NCPU = 1

# make some random data
SEED = 42                    # set a random seed to control
NOBJ = 5                     # number of point sources to simulate
NIMG = 5                     # number of random images
MRANGE = (17., 20.)          # range of magnitude (in FILTER)
APERRAD = 0.3                # radius of the segmentation
PSFSIG = 1.5                 # the PSF sigma
SEDPATH = 'cdbs'


def make_scene():
    np.random.seed(SEED)

    # create something to store the SEDs that will be used for the simulation
    hdul = fits.HDUList()

    # make a direct and segmentation image
    img = np.zeros(SIZE, dtype=float)
    seg = np.zeros_like(img, dtype=int)
    xim, yim = np.meshgrid(np.arange(SIZE[0]), np.arange(SIZE[1]))

    # create a morphological function for the sources
    func = models.Gaussian2D(x_stddev=PSFSIG, y_stddev=PSFSIG)

    # read the filter curve
    band = su.photometry.Throughput.from_keys(TELESCOPE, INSTRUMENT, FILTER)

    # make dir to stage the SEDs from CDBS
    os.makedirs(SEDPATH, exist_ok=True)

    # radius for the segmaps
    radius = APERRAD / PIXSCL

    # process each
    for segid in range(1, NOBJ + 1):

        # make random positions/magnitude for each source
        x = np.random.uniform(low=0, high=SIZE[0])
        y = np.random.uniform(low=0, high=SIZE[0])
        m = np.random.uniform(low=MRANGE[0], high=MRANGE[1])

        # get a spectral file
        idx = np.random.randint(low=1, high=131)
        basename = f'pickles_uk_{idx}.fits'
        localfile = os.path.join(SEDPATH, basename)
        if not os.path.exists(localfile):
            localfile = su.photometry.SED.get_from_CDBS('pickles', basename,
                                                        outpath=SEDPATH)
            if localfile is None:
                warnings.warn("Did not download file.  May cause problems.",
                              RuntimeWarning)

        hdr = fits.Header()
        hdr['EXTNAME'] = (str(segid), 'segmentation ID')
        hdr['EXTVER'] = (0, 'spectral region')
        hdr['FILENAME'] = (localfile, 'filename')
        hdul.append(fits.TableHDU(header=hdr))

        # problem here.  Astropy uses "amplitude" as the overall
        # normalization of the spatial profile, even though what we want
        # is the total flux to be specified.  Of course for some analytic
        # profiles, this amplitude can be determined, however I want
        # to leave the profile flexible, so will determine the conversion
        # between amplitude and total by just scaling via a "delta"
        # image (called dimg)

        # evaluate the point source
        func.x_mean = x
        func.y_mean = y
        func.amplitude = 1.
        dimg = func(xim, yim)

        # do the normalization
        ftot = 10.**(-0.4 * (m - band.zeropoint))

        # scale the delta image and add into the image
        img += (ftot * (dimg / np.sum(dimg)))

        # compute an aperture for the segmap
        rim = np.hypot(xim - x, yim - y)
        seg = np.where(rim < radius, segid, seg)

    # make a WCS
    w = WCS(naxis=2)
    w.wcs.crpix = (SIZE[0] / 2., SIZE[1] / 2.)
    w.wcs.crval = (RA, DEC)
    w.wcs.ctype = ('RA---TAN', 'DEC--TAN')
    w.wcs.cd = ((-PIXSCL / 3600., 0.), (0., PIXSCL / 3600.))

    # make a header
    h = w.to_header()
    h['TELESCOP'] = TELESCOPE
    h['INSTRUME'] = INSTRUMENT
    h['FILTER'] = FILTER
    h['ZERO'] = band.zeropoint

    # save the images
    fits.writeto(f'{ROOT}_seg.fits', seg, header=h, overwrite=True)
    fits.writeto(f'{ROOT}_sci.fits', img, header=h, overwrite=True)
    hdul.writeto(f'{ROOT}_seds.fits', overwrite=True)


def getPrefix(ins):
    if ins == 'WFC3':
        return 'i'
    elif ins == 'ACS':
        return 'j'
    else:
        raise NotImplementedError("Invalid instrument")


def simulate_grisms(tabulate=False):
    ins = getPrefix(INSTRUMENT)

    # make a file that describes the WCS of the grism images to be made
    # nota bene: these columns can be in any order
    wcsfile = f'{ROOT}_wcs.csv'
    with open(wcsfile, 'w') as fp:
        print('dataset,ra,dec,orientat,telescope,instrument,exptime,disperser,blocking', file=fp)
        for i in range(1, NIMG + 1):
            dataset = f'{ins}{i:07}q'
            orientat = i * (180. / NIMG)
            print(dataset, RA, DEC, orientat, TELESCOPE, INSTRUMENT + DETECTOR,
                  EXPTIME, DISPERSER, BLOCKING, sep=',', file=fp)

    # load the grism images
    data = su.wfss.WFSSCollection.from_wcsfile(wcsfile)

    # load the sources
    sources = su.sources.SourceCollection(f'{ROOT}_seg.fits',
                                          f'{ROOT}_sci.fits',
                                          local_back=False,
                                          sedfile=f'{ROOT}_seds.fits')

    # project the sources onto the grism images
    if tabulate:
        tab = su.modules.Tabulate(ncpu=NCPU)
        _ = tab(data, sources)

    # use the projection tables and the SEDs to simulate grism images
    sim = su.modules.Simulate(ncpu=NCPU)
    imgfiles = sim(data, sources)  # noqa: F841


def run_all():

    make_scene()

    simulate_grisms(tabulate=True)


if __name__ == '__main__':
    run_all()
