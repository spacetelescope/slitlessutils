from astropy.io import fits
from astropy.wcs import WCS
from astroquery.mast import Observations
from drizzlepac import astrodrizzle

import matplotlib.pyplot as plt
import numpy as np
import os
import shutil

import slitlessutils as su

'''
  1) download data
     a) copy files to working directory
  2) preprocess grism images
     a) flag cosmic rays
     b) subtract master sky
     c) sync WCS
  3) preprocess direct images
     a) run astrodrizzle
     b) make segmentation map
  4) extract single
'''


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
RA = 264.1019010      # in deg
DEC = -32.9087315     # in deg
RAD = 0.5             # in arcsec
SUFFIX, DRZSUF = 'flc', 'drc'
SCALE = 0.05          # driz image pix scale


def download():
    obs_ids = tuple(d.lower() for k, v in DATASETS.items() for d in v)
    obstab = Observations.query_criteria(obs_id=obs_ids, obs_collection=TELESCOPE)

    # some settings for downloading
    kwargs = {'productSubGroupDescription': [SUFFIX.upper()],
              'extension': 'fits'}

    for row in obstab:
        obsid = str(row['obsid'])
        downloads = Observations.download_products(obsid, **kwargs)
        for local in downloads['Local Path']:
            local = str(local)
            f = os.path.basename(local)
            if f.startswith(obs_ids):
                shutil.copy2(local, '.')


def preprocess_grism():

    # create a background subtraction object
    back = su.core.preprocess.background.Background()

    for imgdset, grismdset in zip(DATASETS[FILTER], DATASETS[GRATING]):
        grismfile = f'{grismdset}_{SUFFIX}.fits'
        imgfile = f'{imgdset}_{SUFFIX}.fits'

        # flag CRs by Laplace filtering
        su.core.preprocess.crrej.laplace(grismfile, inplace=True)

        # subtract background via master-sky
        back.master(grismfile)

        # update WCS to match Gaia
        su.core.preprocess.astrometry.upgrade_wcs(imgfile, grismfile,
                                                  inplace=True)


def preprocess_direct():
    files = []
    for imgdset in DATASETS[FILTER]:
        imgfile = f'{imgdset}_{SUFFIX}.fits'
        # su.core.preprocess.crrej.laplace(imgfile,inplace=True)
        files.append(imgfile)

    # mosaic data via astrodrizzle
    astrodrizzle.AstroDrizzle(files, output=ROOT, build=False,
                              static=False, skysub=True, driz_separate=False,
                              median=False, blot=False, driz_cr=False,
                              driz_combine=True, final_wcs=True,
                              final_rot=0., final_scale=SCALE,
                              final_pixfrac=1.0,
                              overwrite=True, final_fillval=0.0)

    # AGH gotta remove second extensions
    # Must use memmap=False to force close all handles and allow file overwrite
    with fits.open(f'{ROOT}_{DRZSUF}_sci.fits', memmap=False) as hdulist:
        img = hdulist['PRIMARY'].data
        hdr = hdulist['PRIMARY'].header

    wcs = WCS(hdr)
    x, y = wcs.all_world2pix(RA, DEC, 0)

    xx, yy = np.meshgrid(np.arange(hdr['NAXIS1']),
                         np.arange(hdr['NAXIS2']))
    rr = np.hypot(xx-x, yy-y)
    seg = rr < (RAD/SCALE)

    # add some things for SU
    hdr['TELESCOP'] = TELESCOPE
    hdr['INSTRUME'] = INSTRUMENT
    hdr['FILTER'] = FILTER

    # write the files to disk
    fits.writeto(f'{ROOT}_{DRZSUF}_sci.fits', img, hdr, overwrite=True)
    fits.writeto(f'{ROOT}_{DRZSUF}_seg.fits', seg.astype(int), hdr, overwrite=True)


def extract_single():

    # load data into SU
    files = [f'{f}_{SUFFIX}.fits' for f in DATASETS[GRATING]]
    data = su.wfss.WFSSCollection.from_list(files)

    # load the sources into SU
    sources = su.sources.SourceCollection(f'{ROOT}_{DRZSUF}_seg.fits',
                                          f'{ROOT}_{DRZSUF}_sci.fits',
                                          zeropoint=ZEROPOINT)

    # project the sources onto the grism images
    tab = su.modules.Tabulate(ncpu=1)
    pdtfiles = tab(data, sources)  # noqa: F841

    # run the single-orient extraction
    ext = su.modules.Single('+1', mskorders=None, root=ROOT)
    res = ext(data, sources)  # noqa: F841


def plot():

    # l,f=np.loadtxt('wr96_hres.dat',unpack=True,usecols=(0,1))
    # ff=gaussian_filter1d(f,32)

    dat, hdr = fits.getdata(f'{ROOT}_x1d.fits', header=True)
    dat['flam'] *= 1e-17

    plt.plot(dat['lamb'], dat['flam']/1e-13, label=GRATING)
    # plt.plot(l,f/1e-13,label='Larsen et al. (high-res)')
    # plt.plot(l,ff/1e-13,label='Larsen et al. (smoothed)')
    plt.ylabel(r'$f_{\lambda}$ ($10^{-13}$ erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$)')
    plt.xlabel(r'$\lambda$ ($\mathrm{\AA}$)')
    plt.legend()
    plt.ylim(0., 1.4)

    plt.xlim(5500, 10000)

    plt.savefig('wr96.pdf')
    plt.show()


def run_all(plot=True):
    su.start_logging()
    download()
    preprocess_grism()
    preprocess_direct()
    extract_single()
    if plot:
        plot()


if __name__ == '__main__':

    run_all()
