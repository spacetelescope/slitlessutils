import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from scipy.special import gammaincinv

import slitlessutils as su

from .parameters import DEC, NPIX, PIXSCL, RA

# import matplotlib as mpl
# import matplotlib.pyplot as plt
# from matplotlib.patches import Ellipse


ROOT = 'ring'
DATASETS = ('a', 'b', 'c')      # will simulate 3 datasets
TELESCOPE = 'HST'
INSTRUMENT = 'ACS'
DETECTOR = 'WFC'
DISPERSER = 'G800L'
FILTER = 'F775W'              # for the direct image
BLOCKING = ''
SEGID = -1

SPECCAT = 'bc95'
SPECFILE = 'bc95_f_10E9.fits'
SEDPATH = 'ring_input'


def make_scene():

    # source properties
    mag = 18.         # total magnitude in an aperture (AB)
    emm_fact = 1e-2   # multiplier to increase emission line strength
    R = 20            # aperture size along the major axis (in pixels)

    # continuum morphology
    n = 1.           # Sersic index
    Re = 6.2         # effective radius (pix)
    q = 0.47         # axis ratio
    pa = 73.         # position angle (deg)

    # get a continuum spectrum
    contsed = su.photometry.SED.from_CDBS(SPECCAT, SPECFILE)

    # emission line morphology
    emm_rad = 4.1   # centroid of emission line region (pix)
    emm_sig = 0.9    # width of emission line region (pix)

    # emission line spectrum
    clam = 6563.     # central wavelength of Gaussian emission line (A)
    slam = 10.       # sigma of Gaussian Emission line (A)

    # make a grid for the images
    d = NPIX/2

    # establish a Cartesian grid
    y, x = np.mgrid[-d:d+1:1, -d:d+1:1]
    # extent=(-d,d,-d,d)           # define the extent of the image

    # rotate the coordinates
    theta = np.radians(pa)
    sn = np.sin(theta)
    cs = np.cos(theta)
    xx = +x*cs+y*sn
    yy = -x*sn+y*cs

    # transform to elliptical coordinates
    rad = np.hypot(xx, yy/q)

    # make images that describe relative contribution of continuum
    # and emission line

    # Sersic value
    bn = gammaincinv(2.*n, 0.5)

    # normalization of the continuum image
    contnorm = np.exp(-bn*((rad/Re)**(1./n)-1.))
    contnorm /= np.sum(contnorm)

    # normalization of the emission line image
    emmnorm = np.exp(-0.5*((rad-emm_rad)/emm_sig)**2)
    emmnorm /= np.sum(emmnorm)

    # read the filter
    band = su.photometry.Throughput.from_keys(TELESCOPE, INSTRUMENT, FILTER)

    # emission line spectrum (assumed as a Gaussian)
    emmflam = np.exp(-0.5*((contsed.lamb-clam)/slam)**2)
    emmflam /= np.sqrt(2*np.pi*slam*slam)
    emmsed = su.photometry.SED(contsed.lamb, emmflam)

    # compute bandpass averaged integrals over continuum and
    # emission line spectra
    fcnt = su.photometry.avefnu(contsed, band)
    femm = su.photometry.avefnu(emmsed, band)

    # make a total image as a sum ofer components, weighted by their
    # bandpass-weighted flux.  Note: the `emm_fact` is the relative
    # scaling between continuum and emission line
    gal = contnorm*fcnt+emmnorm*femm*emm_fact

    # now only take the points inside an aperture
    g = np.where(rad < R)
    # b = np.where(rad >= R)
    fint = np.sum(gal[g])   # instrumental flux (PRE-SCALE)
    # fout=np.sum(gal[b])
    # ftot=fint+fout

    scl = 10**(-0.4*(mag+48.6))/fint    # scale factor to shift fluxes

    # apply scale factor to images and shift units to e-
    # units=(2.998e10/band.photplam)*(1e8/band.photplam)/band.photflam

    cntimg = contnorm*fcnt*scl/band.photfnu
    emmimg = emmnorm*femm*scl/band.photfnu

    # the final galaxy image is a sum over the spectral components
    galimg = cntimg+emmimg

    # make a segimage
    segimg = np.zeros_like(galimg, dtype=int)
    segimg[g] = SEGID         # make the SEGID negative, for pix-by-pix decomp

    # make an SED file
    hdul = fits.HDUList()
    for regid, (y, x) in enumerate(zip(*g)):
        sed = contsed*fcnt*contnorm[y, x]+emmsed*femm*emmnorm[y, x]*emm_fact
        hdu = sed.as_HDU(extname=str(SEGID), extver=regid,
                         segid=SEGID, regid=regid, x=x, y=y)

        hdul.append(hdu)
    hdul.writeto(f'{ROOT}_seds.fits', overwrite=True)

    # make a header
    w = WCS(naxis=2)
    w.wcs.crpix = [d, d]
    w.wcs.crval = [RA, DEC]
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    w.wcs.cd = [[-PIXSCL/3600., 0.], [0., PIXSCL/3600.]]
    h = w.to_header()

    # add some info to the file
    h['TELESCOP'] = TELESCOPE
    h['INSTRUME'] = INSTRUMENT
    h['FILTER'] = FILTER
    h['ZERO'] = band.zeropoint

    # write the images to disk
    fits.writeto(f'{ROOT}_seg.fits', segimg, header=h, overwrite=True)
    fits.writeto(f'{ROOT}_sci.fits', galimg, header=h, overwrite=True)


def simulate_grisms():

    # write a "WCS" file to disk that contains properties of
    # the images to emulate an observers setup
    with open(f'{ROOT}_wcs.csv', 'w') as fp:
        print('dataset,ra,dec,orientat,telescope,instrument,disperser,blocking', file=fp)

        n = len(DATASETS)

        # put each WFSS exposure in the canon
        for i, dset in enumerate(DATASETS):
            orientat = i*(180./n)
            line = ','.join((dset+'_'+ROOT, str(RA), str(DEC), str(orientat),
                             TELESCOPE, INSTRUMENT+DETECTOR, DISPERSER, BLOCKING))
            print(line, file=fp)

    # load the grism images
    data = su.wfss.WFSSCollection.from_wcsfile(f'{ROOT}_wcs.csv')

    # load the sources
    sources = su.sources.SourceCollection(f'{ROOT}_seg.fits', f'{ROOT}_sci.fits',
                                          sedfile=f'{ROOT}_seds.fits')

    # save the normalized SEDs
    sources.write_seds(SEDPATH)

    # project the sources onto the grism images
    # tab = su.modules.Tabulate()
    # pdtfiles = tab(data, sources)

    # use the projection tables and the SEDs to simulate grism images
    sim = su.modules.Simulate(ncpu=1)
    imgfiles = sim(data, sources)

    return imgfiles


def run_all():
    """
    Method to run all of the substeps in a single call
    (Order copied from starfield example)

    Notes
    -----
    1) See individual functions for their descriptions

    """

    make_scene()
    simulate_grisms()
