import os
import shutil

from astropy.io import fits
from astropy.wcs import WCS
from astroquery.mast import Observations
from drizzlepac import astrodrizzle
import matplotlib.pyplot as plt
import numpy as np

import slitlessutils as su


"""
Description
-----------

This file demonstrates the analysis of HST WFC3/UVIS grism spectroscopy using
the G280 observations of the supernova PTF12DAM.  It will compare the
spectrum extracted with slitlessutils to the same measurements from HSTaXe
that is stored in the reference file database.

Author
------
R. Ryan (STScI)

Date
----
Nov 6, 2025
"""


# the observations
TELESCOPE = 'HST'
INSTRUMENT = 'WFC3'
GRATING = 'G280'
FILTER = 'F200LP'
ZEROPOINT = 25.655    # a dummy variable

# properties of source
RA = 216.1925314
DEC = 46.2301076
RAD = 0.4     # in arcsec

# properties of the drizzle products
ROOT = 'PTF12DAM'
SUFFIX, DRZSUF = 'flc', 'drc'
SCALE = 0.04          # driz image pix scale

# datasets to download
DATASETS = {GRATING: ["ibr502010"],
            FILTER: ["ibr502i6q"]}


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
    # build background subtraction
    back = su.core.preprocess.background.Background()  # noqa: F841

    # process each file
    grismfiles = []
    with open(f'{GRATING}.lst', 'r') as fp:
        for line in fp:
            grismfile = line.strip()

            # roll back the WCS to match the Direct imaging
            su.core.preprocess.astrometry.downgrade_wcs(grismfile, key='A',
                                                        inplace=True)

            grismfiles.append(grismfile)

    # mask the CRs
    su.core.preprocess.crrej.drizzle(grismfiles, grouping=None)


def preprocess_direct():

    imagefiles = []

    with open(f'{FILTER}.lst', 'r') as fp:
        for line in fp:
            imagefile = line.strip()

            # roll back the WCS to match the grism
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

    make_segmap()

    print("DONE\n\n\n\n\n\n\n")


def make_segmap():

    # gotta remove second extensions
    # Must use memmap=False to force close all handles and allow file overwrite
    with fits.open(f'{ROOT}_{DRZSUF}_sci.fits', memmap=False) as hdulist:
        img = hdulist['PRIMARY'].data
        hdr = hdulist['PRIMARY'].header

    wcs = WCS(hdr)
    x, y = wcs.all_world2pix(RA, DEC, 0)

    xx, yy = np.meshgrid(np.arange(hdr['NAXIS1']),
                         np.arange(hdr['NAXIS2']))
    rr = np.hypot(xx - x, yy - y)
    seg = rr < (RAD / SCALE)

    # add some things for SU
    hdr['TELESCOP'] = TELESCOPE
    hdr['INSTRUME'] = INSTRUMENT
    hdr['FILTER'] = FILTER

    # write the files to disk
    fits.writeto(f'{ROOT}_{DRZSUF}_sci.fits', img, hdr, overwrite=True)
    fits.writeto(
        f'{ROOT}_{DRZSUF}_seg.fits',
        seg.astype(int),
        hdr,
        overwrite=True)


def extract_single(regions=False, tabulate=True):

    # load data into SU
    # files = [f'{f}_{SUFFIX}.fits' for f in DATASETS[GRATING]]
    data = su.wfss.WFSSCollection.from_file(f'{GRATING}.lst')

    # load the sources into SU
    sources = su.sources.SourceCollection(f'{ROOT}_{DRZSUF}_seg.fits',
                                          f'{ROOT}_{DRZSUF}_sci.fits',
                                          local_back=False,
                                          zeropoint=ZEROPOINT)

    if regions:
        reg = su.modules.Region()
        regfiles = reg(data, sources)  # noqa: F841

    # project the sources onto the grism images
    if tabulate:
        tab = su.modules.Tabulate(ncpu=1)
        pdtfiles = tab(data, sources)  # noqa: F841

    # run the single-orient extraction
    ext = su.modules.Single('+1', mskorders=None, root=ROOT, ncpu=1)
    res = ext(data, sources, profile='uniform')  # noqa: F841


def plot_spectra():
    cfg = su.config.Config()

    # read the SU spectrum
    dat, hdr = fits.getdata(f'{ROOT}_x1d.fits', header=True)

    # fetch the HSTaXe spectrum
    subpath = os.path.join('instruments', 'WFC3UVIS')
    axefile = cfg.get_reffile('ibr502i9q_flc_2.SPC.fits', subpath)

    # read the HSTaXe spectrum
    axe = fits.getdata(axefile, ext=3)
    axe = axe[np.argsort(axe['LAMBDA'])]
    axe['FLUX'] /= cfg.fluxscale
    axe['FERROR'] /= cfg.fluxscale

    # make the plot
    plt.plot(dat['lamb'], dat['flam'], label='slitlessutils')
    plt.plot(axe['LAMBDA'], axe['FLUX'], label='aXe')

    plt.legend()

    plt.xlim(1900, 7700)
    plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')

    plt.ylim(20, 200)
    unit = r'$10^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$)'
    plt.ylabel(r'$f_{\lambda}$ (' + unit + ')')

    plt.tight_layout()
    plt.savefig('wfc3uvis_ptf12dam.pdf')
    plt.show()


def run_all(plot=True):
    # Step 1: Fetch the data
    download()

    # Step 2: Preprocess the grism images
    preprocess_grism()

    # Step 3: Preprocess the direct images
    preprocess_direct()

    # Step 4 (optional):  List the astrometry to verify that everything
    # has the *SAME* WCSNAME
    su.core.preprocess.astrometry.list_wcs('i*flc.fits')

    # Step 5.  Extract the 1d spectra
    extract_single()

    # Step 6 (optional):  Plot the 1d spectra
    if plot:
        plot_spectra()


if __name__ == '__main__':

    run_all()
