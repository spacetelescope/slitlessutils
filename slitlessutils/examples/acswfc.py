import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.modeling import fitting, models
from astropy.wcs import WCS
from astroquery.mast import Observations
from drizzlepac import astrodrizzle
from scipy.ndimage import gaussian_filter1d

import slitlessutils as su

"""
Description
-----------

This file demonstrates the analysis of HST ACS/WFC grism spectroscopy using
the G800L observations of the standard star WR96.  It will compare the
spectrum extracted with slitlessutils to a reference spectrum from
Pasquali et al. (2002) that is stored in the reference file database.

The grey bar shows the region used to normalize the spectra.

Author
------
R. Ryan (STScI)

Date
----
Nov 6, 2025
"""


# the observations
TELESCOPE = 'HST'
INSTRUMENT = 'ACS'
GRATING = 'G800L'
FILTER = 'F775W'
ZEROPOINT = 25.655

# datasets to process
DATASETS = {GRATING: ["jdql01jpq", "jdql01jxq"],
            FILTER: ["jdql01jnq", "jdql01jvq"]}

# output root for astrodrizzle
ROOT = 'WRAY-15-1736'

# position of source
RA = 264.1018744
DEC = -32.9087434
RAD = 0.5             # in arcsec
SUFFIX, DRZSUF = 'flc', 'drc'
SCALE = 0.05          # driz image pix scale


def download():
    # some settings for downloading
    kwargs = {'productSubGroupDescription': [SUFFIX.upper()],
              'extension': 'fits'}

    for band in DATASETS.keys():
        with open(f'{band}.lst', 'w') as fp:

            for dataset in DATASETS[band]:
                ipppss = dataset[:6]
                obstab = Observations.query_criteria(
                    obs_id=dataset, obs_collection=TELESCOPE)

                for row in obstab:
                    obsid = str(row['obsid'])
                    downloads = Observations.download_products(obsid, **kwargs)
                    for local in downloads['Local Path']:
                        filename = os.path.basename(local)
                        if filename.startswith(ipppss):

                            shutil.copy2(local, '.')
                            print(filename, file=fp)


def preprocess_grism():

    # create a background subtraction object
    back = su.core.preprocess.background.Background()

    # process each image
    grismfiles = []
    with open(f'{GRATING}.lst') as fp:
        for line in fp:
            grismfile = line.strip()

            # subtract background via master-sky
            back.master(grismfile, inplace=True)

            # downgrade the WCS to be consistent
            su.core.preprocess.astrometry.downgrade_wcs(grismfile, key='A',
                                                        inplace=True)

            grismfiles.append(grismfile)

    # flag the cosmic rays
    su.core.preprocess.crrej.drizzle(grismfiles, grouping=None)


def preprocess_direct():

    imagefiles = []
    with open(f'{FILTER}.lst') as fp:
        for line in fp:
            imagefile = line.strip()

            # make the WCS do the older version
            su.core.preprocess.astrometry.downgrade_wcs(imagefile, key='A',
                                                        inplace=True)

            imagefiles.append(imagefile)

    # mosaic data via astrodrizzle
    astrodrizzle.AstroDrizzle(imagefiles, output=ROOT, build=False,
                              static=False, skysub=True, driz_separate=False,
                              median=False, blot=False, driz_cr=False,
                              driz_combine=True, final_wcs=True,
                              final_rot=0., final_scale=SCALE,
                              final_pixfrac=1.0,
                              overwrite=True, final_fillval=0.0)

    # gotta remove second extensions
    # Must use memmap=False to force close all handles and allow file overwrite
    with fits.open(f'{ROOT}_{DRZSUF}_sci.fits', memmap=False) as hdulist:
        img = hdulist['PRIMARY'].data
        hdr = hdulist['PRIMARY'].header

    wcs = WCS(hdr)
    x, y = wcs.all_world2pix(RA, DEC, 0)
    xim, yim = np.meshgrid(np.arange(hdr['NAXIS1']),
                           np.arange(hdr['NAXIS2']))

    sz = 10
    xmin = int(x) - sz
    xmax = int(x) + sz
    ymin = int(y) - sz
    ymax = int(y) + sz
    pix = (slice(ymin, ymax), slice(xmin, xmax))

    mod = models.Gaussian2D(x_mean=sz, y_mean=sz, amplitude=np.amax(img[pix]))

    fitter = fitting.LevMarLSQFitter()
    yy, xx = np.indices(img[pix].shape, dtype=float)

    res = fitter(mod, xx, yy, img[pix])
    x = res.x_mean.value + xmin
    y = res.y_mean.value + ymin

    rr = np.hypot(xim - x, yim - y)
    seg = rr < (RAD / SCALE)
    seg = seg.astype(int)

    # add some things for SU
    hdr['TELESCOP'] = TELESCOPE
    hdr['INSTRUME'] = INSTRUMENT
    hdr['FILTER'] = FILTER

    # write the files to disk
    fits.writeto(f'{ROOT}_{DRZSUF}_sci.fits', img, hdr, overwrite=True)
    fits.writeto(f'{ROOT}_{DRZSUF}_seg.fits', seg, hdr, overwrite=True)


def extract_single(tabulate=True):

    # load data into SU
    # files = [f'{f}_{SUFFIX}.fits' for f in DATASETS[GRATING]]
    # data = su.wfss.WFSSCollection.from_list(files)
    data = su.wfss.WFSSCollection.from_file(f'{GRATING}.lst')

    # load the sources into SU
    sources = su.sources.SourceCollection(f'{ROOT}_{DRZSUF}_seg.fits',
                                          f'{ROOT}_{DRZSUF}_sci.fits',
                                          zeropoint=ZEROPOINT)

    # project the sources onto the grism images
    if tabulate:
        tab = su.modules.Tabulate(ncpu=1, orders=('+1',), remake=True)
        pdtfiles = tab(data, sources)  # noqa: F841

    # run the single-orient extraction
    ext = su.modules.Single('+1', mskorders=None, root=ROOT, ncpu=1)
    res = ext(data, sources, profile='uniform')  # noqa: F841


def plot_spectra():

    cfg = su.config.Config()

    # get the Larsen et al. reference spectrum
    subpath = os.path.join('instruments', 'ACSSBC')
    reffile = cfg.get_reffile('wr96_hres.dat', subpath)

    # load and smooth the Larsen spectrum
    l, f = np.loadtxt(reffile, unpack=True, usecols=(0, 1))
    f /= 1e-13

    # smooth the A. Pasquali reference spectrum (kindly provided S. Larsen).
    # the Scale factor is from comparing the notional ACS dispersion
    # (40A/pix) to the spectrum quoted from Pasquali+ 2002 which has
    # 1.26 A/pix. Therefore smoothing factor is 40./1.26 = 31.7
    ff = gaussian_filter1d(f, 31.7)

    # load the data and change the units
    dat, hdr = fits.getdata(f'{ROOT}_x1d.fits', header=True)
    dat['flam'] *= cfg.fluxscale / 1e-13
    dat['func'] *= cfg.fluxscale / 1e-13

    # get a good range of points to compute a (variance-weighted) scale factor
    lmin = 7350.
    lmax = 8120.
    g = np.where((lmin <= dat['lamb']) & (dat['lamb'] <= lmax))[0]
    ff2 = np.interp(dat['lamb'], l, ff)
    obssn = dat['flam'][g] / dat['func'][g]
    calsn = ff2[g] / dat['func'][g]
    den = np.nansum(obssn * obssn)
    num = np.nansum(calsn * obssn)
    scl = num / den

    # we don't expect the scales to match, since the apertures of the
    # Pasquali spectrum match ours
    print(f'Scaling factor: {scl}')

    # plot the SU spectrum
    plt.axvspan(lmin, lmax, color='lightgrey')
    plt.plot(dat['lamb'], dat['flam'], label='slitlessutils')
    plt.plot(l, ff / scl, label='Pasquali et al. (2002) $-$ smoothed')

    # uncomment this to see the hi-res file.
    # plt.plot(l, f, label='Pasquali et al. (high-res)')

    # label the axes
    unit = r'$10^{-13}$ erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$)'
    plt.ylabel(r'$f_{\lambda}$ (' + unit + ')')
    plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')

    # put on legend and change limits
    plt.legend(loc='upper left')
    plt.ylim(0.0, 1.2)
    plt.xlim(5500, 10000)
    plt.tight_layout()

    # write the file to disk
    plt.savefig('acswfc_wr96.pdf')
    plt.show()


def run_all(plot=True):

    # step 1.  Fetch the data
    download()

    # step 2.  Process the grism images
    preprocess_grism()

    # step 3.  process the direct images
    preprocess_direct()

    # step 4 (optional).  List the astrometry to verify that everything
    # has the *SAME* WCSNAME
    su.core.preprocess.astrometry.list_wcs('j*flc.fits')

    # step 5.  Extract the 1d spectra
    extract_single(tabulate=True)

    # step 6 (optional).  Plot the 1d spectra
    if plot:
        plot_spectra()


if __name__ == '__main__':

    run_all()
