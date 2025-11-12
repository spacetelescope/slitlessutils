import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import sigma_clipped_stats

from .....config import Config
from ....utilities import indices
from .spectraltable import SpectralTable


def boxcar(source, det, sci, unc, dqa, flatfield, order, odt,
           width=17, summation='cartesian', plot=False):
    '''
    Worker function to perform boxcar or "aperture" spectroscopy

    Parameters
    ----------
    source : `slitlessutils.core.sources.Source`
       The source for which to do the extraction

    det : `slitlessutils.core.wfss.data.WFSSDetector()`
       The data for which to do the extraction

    sci : `np.ndarray`
       The science image to do the extraction on.

    unc : `np.ndarray`
       The uncertainty image to do the extraction on.

    dqa : `np.ndarray`
       The data-quality image to do the extraction on.

    flatfield : `slitlessutils.core.wfss.config.Flatfield`
       The flat-field object to use.

    order : `slitlessutils.core.wfss.config.Order`
       The spectral order to use.

    odt : `slitlessutils.core.tables.ODT`
       The object-dispersion table (ODT).

    width : int or float
       The extraction width in pixels.  Default is 15 pix

    summation : str
       The mode of the pixel summation.  Default is 'cartesian'.

    plot : bool
       Flag to make a diagnostic plot.  Default is False

    Returns
    -------
    spectrum : `SpectralTable`
       The measured spectrum with ancillary data.

    '''

    ny, nx = sci.shape

    fluxscale = Config().fluxscale

    # distill out properties
    xg = odt.get('x') + 1
    yg = odt.get('y') + 1
    vg = odt.get('val')
    wav = odt.wavelengths()

    # make an output data structure to fill
    spec = SpectralTable(odt.segid)

    # half width of the aperture
    half = (width - 1) / 2

    # get the center
    x0, y0 = det.ad2xy(*source.adc)

    # in fixed stripes (ie. Cartesian)
    if summation == 'cartesian':
        # find unique pixels
        xi = indices.reverse(xg)
        # yi = indices.reverse(yg)

        ytrace = np.zeros(len(xi))
        xtrace = np.zeros(len(xi))
        for i, (x, idx) in enumerate(xi.items()):

            # collect the data for this slice
            # xu = xg[idx]
            yu = yg[idx]
            vu = vg[idx]
            wu = wav[idx]

            # compute some things for this slice
            wave = np.average(wu, weights=vu)
            ytrace[i] = np.average(yu, weights=vu)
            xtrace[i] = x

            # get the fractional part of the trace to use for
            # fractional pixel coverage
            ycent, yfrac = np.divmod(ytrace[i], 1)

            # get the range
            ymin = max(int(ycent - half), 0)
            ymax = min(int(ycent + half) + 1, ny - 1)

            # create dummy arrays to vectorize calculations
            ys = np.arange(ymin, ymax + 1, 1, dtype=int)
            xs = np.full_like(ys, x, dtype=int)

            # make extraction aperture (always length = 1+width)
            aper = np.ones_like(ys, dtype=float)
            aper[0] = 1 - yfrac
            aper[-1] = yfrac
            npix = np.sum(aper)

            # compute calibrations
            disp = order.dispersion(x0, y0, wavelength=wave)
            sens = order.sensitivity(wave) * fluxscale
            flat = flatfield(xs, ys, wave)
            area = det.relative_pixelarea(xs, ys)
            calib = flat * area * sens * disp

            # apply the calibrations
            maxcal = np.amax(calib)
            if maxcal > 1e-15:
                # the calibrated images
                calsci = (sci[ys, xs] / calib)
                calvar = (unc[ys, xs] / calib)**2

                # sum the data with calibrations applied
                flam = np.sum(calsci * aper)
                fvar = np.sum(calvar * aper)
            else:
                flam = np.inf
                fvar = np.inf

            # sum the data
            flux = np.nansum(sci[ys, xs] * aper)

            # save the results
            spec.append(wave, disp, flam, np.sqrt(fvar), flux, 0.0, npix)

        # make a diagnostic plot?
        if plot:
            bbx = slice(np.amin(xg), np.amax(xg))
            bby = slice(np.amin(yg), np.amax(yg))

            vv, yy, xx = indices.decimate(vg, yg, xg)
            mod = np.zeros_like(sci)
            mod[yy, xx] = vv

            fig, ax = plt.subplots(2, 1, sharex=True, sharey=True,
                                   figsize=(10, 5), layout='constrained')
            ax[0].imshow(mod, vmin=0, vmax=1)
            ax[0].errorbar(xtrace, ytrace, yerr=np.full_like(ytrace, width / 2),
                           color='red')
            a, m, s = sigma_clipped_stats(sci[bby, bbx])

            ax[1].imshow(sci, vmin=m - 3 * s, vmax=m + 3 * s)
            ax[1].errorbar(xtrace, ytrace, yerr=np.full_like(ytrace, width / 2),
                           color='red')

            ax[0].set_xlim(np.amin(xtrace), np.amax(xtrace))
            ax[0].set_ylim(np.amin(ytrace) - 20, np.amax(ytrace) + 20)

            plt.show()

    elif summation == 'perpendicular':
        raise NotImplementedError("Not yet supported.")

    return spec
