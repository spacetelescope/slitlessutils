import os
import shutil

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astroquery.mast import Observations
from drizzlepac import astrodrizzle

import slitlessutils as su

"""
Description
-----------

This file demonstrates the analysis of HST WFC3/IR grism spectroscopy using
the G102 observations of the standard star GD153.  It will compare the
spectrum extracted with slitlessutils to the reference spectrum, taken
from CALSPEC but redistributed with slitlessutils.

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
GRATING = 'G102'
FILTER = 'F105W'
ZEROPOINT = 26.9

# datasets to process
DATASETS = {FILTER: ['id2q01flq'],
            GRATING: ['id2q01fmq', 'id2q01fpq', 'id2q01fsq']}

# output root for astrodrizzle
ROOT = 'GD153'

# position of source
RA = 194.2595160
DEC = 22.0304329
RAD = 2.0             # in arcsec
# RAD = 0.5
SUFFIX, DRZSUF = 'flt', 'drz'
SCALE = 0.12          # driz image pix scale


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
    with open(f'{GRATING}.lst') as fp:
        for line in fp:
            grismfile = line.strip()

            # subtract the background
            back.master(grismfile, inplace=True)

            # update the astrometry
            su.core.preprocess.astrometry.downgrade_wcs(grismfile, key='A',
                                                        inplace=True)


def preprocess_direct():
    imagefiles = []
    with open(f'{FILTER}.lst') as fp:
        for line in fp:
            imagefile = line.strip()

            # reset the astrometry
            su.core.preprocess.astrometry.downgrade_wcs(
                imagefile, key='A', inplace=True)
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


def extract_single(tabulate=False):

    # load data into SU
    data = su.wfss.WFSSCollection.from_file(f'{GRATING}.lst')

    # load the sources into SU
    sources = su.sources.SourceCollection(f'{ROOT}_{DRZSUF}_seg.fits',
                                          f'{ROOT}_{DRZSUF}_sci.fits',
                                          local_back=False,
                                          zeropoint=ZEROPOINT)

    reg = su.modules.Region()
    regfiles = reg(data, sources)  # noqa: F841

    # project the sources onto the grism images
    if tabulate:
        tab = su.modules.Tabulate(ncpu=1)
        pdtfiles = tab(data, sources)  # noqa: F841

    # run the single-orient extraction
    ext = su.modules.Single('+1', mskorders=None, root=ROOT, ncpu=1)
    res = ext(data, sources, profile='uniform', width=15)  # noqa: F841


def extract_multi():
    raise NotImplementedError("This function is not validated yet.")

    # load data into SU
    data = su.wfss.WFSSCollection.from_files(f'{GRATING}.lst')

    # load the sources into SU
    sources = su.sources.SourceCollection(f'{ROOT}_{DRZSUF}_seg.fits',
                                          f'{ROOT}_{DRZSUF}_sci.fits',
                                          local_back=False,
                                          zeropoint=ZEROPOINT)

    # run the multi-orient extraction
    ext = su.modules.Multi('+1', (-3., 1., 0.1), algorithm='grid')
    res = ext(data, sources, root=ROOT + '_multi')

    return res


def plot_spectra():

    cfg = su.config.Config()

    # get the sensitivity curve just for fun
    subpath = os.path.join('instruments', 'WFC3IR')
    sens_file = cfg.get_reffile('WFC3.IR.G102.1st.sens.2.fits', subpath)
    sens = fits.getdata(sens_file)

    # read the SU spectrum
    dat, hdr = fits.getdata(f'{ROOT}_x1d.fits', header=True)

    # fetch the calspec spectrum for comparison
    csfile = cfg.get_reffile('gd153_mod_012.fits', subpath)
    cs = fits.getdata(csfile)
    cs['FLUX'] /= cfg.fluxscale

    # compute the relative offset
    lrange = (8000, 11550)
    g = np.where((lrange[0] < dat['lamb']) & (dat['lamb'] < lrange[1]))
    fint = np.interp(dat['lamb'][g], cs['WAVELENGTH'], cs['FLUX'])
    num = np.sum((dat['flam'][g] / dat['func'][g]) * (fint / dat['func'][g]))
    den = np.sum((dat['flam'][g] / dat['func'][g])
                 * (dat['flam'][g] / dat['func'][g]))
    scl = num / den
    dif = np.abs(scl - 1) * 100.

    # this should be ~5% and is a known issue with the CALSPEC spectrum
    print(f"Relative offset: {dif:.1f}%")

    # make some plots
    figsize = (6, 6)
    fig, ax = plt.subplots(
        2, 1, sharex=True, layout='constrained', figsize=figsize)

    # plot grey regions showing were to ignore
    color = 'grey'
    alpha = 0.2

    ax[0].axvspan(0, lrange[0], color=color, alpha=alpha)
    ax[0].axvspan(lrange[1], 1e10, color=color, alpha=alpha)
    ax[1].axvspan(0, lrange[0], color=color, alpha=alpha)
    ax[1].axvspan(lrange[1], 1e10, color=color, alpha=alpha)

    # plot the two spectra
    ax[0].errorbar(
        dat['lamb'],
        dat['flam'],
        yerr=dat['func'],
        label='slitlessutils')
    ax[0].plot(cs['WAVELENGTH'], cs['FLUX'], label='CALSPEC')

    # add a legend
    ax[0].legend(loc='upper right')

    # set y axis labels/limits
    ax[0].set_ylim(0, 600)
    unit = r'$10^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}^{-1}$'
    ax[0].set_ylabel(r'$f_{\lambda}$ (' + unit + ')')

    # plot the sensitivity curve in a lower panel to show the wavelength range
    # of confidence
    color = 'red'
    ax[1].plot(
        sens['WAVELENGTH'],
        sens['SENSITIVITY'] / 1e17,
        label='sensitivity',
        color=color)
    ax[1].set_ylabel(
        r'$S_\lambda$ ($10^{17}$ $e^-$/s per erg/cm$^2$/s/$\mathrm{\AA}$)')
    ax[1].set_ylim(0, 0.11)

    # set the x axis labels/limits
    ax[1].set_xlim(7500, 12000)
    ax[1].set_xlabel(r'$\lambda$ ($\mathrm{\AA}$)')

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
    su.core.preprocess.astrometry.list_wcs('i*flt.fits')

    # Step 5.  Extract the 1d spectra
    extract_single(tabulate=True)

    # Step 6 (optional):  Plot the 1d spectra
    if plot:
        plot_spectra()


if __name__ == '__main__':

    run_all()
