from astropy.io import fits
from astropy.wcs import WCS
from astroquery.mast import Observations
from drizzlepac import astrodrizzle
import numpy as np
import pytest


import slitlessutils as su


@pytest.mark.remote_data
@pytest.mark.usefixtures('_jail')
def test_ACS_grism(tmp_path):
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

    obs_ids = tuple(d.lower() for k, v in DATASETS.items() for d in v)
    obstab = Observations.query_criteria(obs_id=obs_ids, obs_collection=TELESCOPE)

    # some settings for downloading
    kwargs = {'productSubGroupDescription': [SUFFIX.upper()],
              'extension': 'fits'}

    for row in obstab:
        obsid = str(row['obsid'])
        Observations.download_products(obsid, **kwargs)

    # create a background subtraction object
    back = su.core.preprocess.background.Background()

    for imgdset, grismdset in zip(DATASETS[FILTER], DATASETS[GRATING]):
        grismfile = f'mastDownload/HST/{grismdset}/{grismdset}_{SUFFIX}.fits'
        imgfile = f'mastDownload/HST/{imgdset}/{imgdset}_{SUFFIX}.fits'

        # flag CRs by Laplace filtering
        su.core.preprocess.crrej.laplace(grismfile, inplace=True)

        # subtract background via master-sky
        back.master(grismfile)

        # update WCS to match Gaia
        su.core.preprocess.astrometry.upgrade_wcs(imgfile, grismfile,
                                                  inplace=True)

    files = []
    for imgdset in DATASETS[FILTER]:
        imgfile = f'{imgdset}_{SUFFIX}.fits'
        # su.core.preprocess.crrej.laplace(imgfile,inplace=True)
        files.append(f"mastDownload/HST/{imgdset}/{imgfile}")

    # mosaic data via astrodrizzle

    astrodrizzle.AstroDrizzle(files, output=f'{ROOT}', build=False,
                              static=False, skysub=True, driz_separate=False,
                              median=False, blot=False, driz_cr=False,
                              driz_combine=True, final_wcs=True,
                              final_rot=0., final_scale=SCALE,
                              final_pixfrac=1.0, preserve=False,
                              overwrite=True, final_fillval=0.0)

    # Have to remove second extensions
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
    fits.writeto(f'{ROOT}_{DRZSUF}_seg.fits', seg.astype(int), hdr, overwrite=True)

    # load data into SU
    files = [f'mastDownload/HST/{f}/{f}_{SUFFIX}.fits' for f in DATASETS[GRATING]]
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

    # load the data and change the units
    cfg = su.config.Config()
    dat, hdr = fits.getdata(f'{ROOT}_x1d.fits', header=True)
    dat['flam'] *= cfg.fluxscale / 1e-13
    dat['func'] *= cfg.fluxscale / 1e-13

    assert dat.shape == (119,)
    assert np.allclose(dat[1], (5620.0, 0.029815594, 0.00016166682, 0.0, 13))
    assert np.allclose(dat[-1], (10300.0, 0.081851676, 0.0007620386, 0.0, 4))
