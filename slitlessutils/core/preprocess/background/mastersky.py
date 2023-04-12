from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import numpy as np
import os

from .cartesian import cartesian
from ....config import Config
from ....logger import LOGGER
from ...utilities import headers
from ...wfss import WFSS


# def get_metadata(filenames):
#    ''' helper function to get metadata from a list of fits filenames '''
#
#
#    visits={}
#
#    for f in filenames:
#        h0=fits.getheader(f,0)
#        dataset=h0.get("ROOTNAME")
#        vis=dataset[4:6]
#        if vis not in visits:
#            visits[vis]=[]
#
#        visits[vis].append(f)
#
#    return visits
#

def mastersky(data, nsig=3.0, epsilon=1e-5, maxiter=10,
              newfile=None, inplace=True):
    """
    Function to fit a *single* master sky image to a WFSS image.

    The master-sky image is a two-dimensional model of the dispersed
    sky background.  Ideally, this image is scaled to match the background
    pixels in a dispersed image, however requires flagging all of the
    pixels that contain any additional signal (such as, but not limited to,
    astrophysical sources from any spectral order, cosmic rays,
    persistence or any type of bad pixels).  The masking is implemented
    via an iterative approach:

    0) initalize the background model as the 3-sigma clipped median
       (:math:`m`) of the pixels not flagged in the data-quality array
       and flag pixels that are

    .. math::
       \\frac{{\\cal I}-m}{{\\cal U}} \\leq n_{sig}

    where :math:`{\\cal I}` and :math:`{\\cal U}` are the science and
    uncertainty images, respectively.

    1) compute optimal scaling coefficient as

    .. math::
       \\alpha = \\frac{F_{OT}}{F_{TT}}

    where

    .. math::
       :nowrap:
       \begin{eqnarray}
         F_{OO} & = & \\sum_{x,y\\in s} \\left(\\frac{{\\cal I}}{{\\cal U}}\\right)^2\\
         F_{OT} & = & \\sum_{x,y\\in s} \\frac{{\\cal I}}{{\\cal U}}\\frac{{\\cal S}}{{\\cal U}}\\
         F_{TT} & = & \\sum_{x,y\\in s} \\left(\\frac{{\\cal S}}{{\\cal U}}\\right)^2\\
       \\end{eqnarray}

    and :math:`s` is the collection of pixels that are neither flagged in
    the data-quality arrays nor above the background level.

    2) redefine the non-background pixels by testing:
    .. math::
       \\frac{{\\cal I}-\\alpha{\\cal S}}{{\\cal U}} \\leq n_{sig}

    3) go to step 1) until either the maximum number of iterations is reached
    or the fractional difference is less than the set tolerance
    (:math:`\\epsilon`):

    .. math::
       |\\alpha_{i-1} - \\alpha_i| \\leq \\epsilon |\\alpha_{i}|

    Parameters
    ----------
    data : str or `wfss.WFSS`
        Either a string that is the full path to a WFSS file or a
        `wfss.WFSS`.

    nsig : float, optional
        The number of sigma above the uncertainty map to identify
        pixels that contain any signals other than only background.
        Default is 3.0

    epsilon : float, optional
        The tolerance for convergence of identifying background vs.
        non-background pixels.  This is given as a fractional change
        of the scaling coeffiecient between successive iterations:
        Default is 1e-5.

    maxiter : int, optional
        Maximum number of iterations before stopping.  Default is 10.

    newfile : str, optional
        The name of the new file.  If this is a valid string, then the
        master-sky-subtracted file will be written to this filename.
        If it is set to `None`, then a name will be created based on the
        dataset of the input file:

        f'{dataset}_msky.fits'

        Default is `None`.

    inplace : bool, optional
        A flag that the outputs should be overwritten to the input file.
        Default is `True`.


    Returns
    -------
    outfile : str
        The name of the master-sky subtracted file.


    Notes
    -----
    This will also write out a source-mask file that is a boolean image
    flagging the pixels that were identified as source-pixels (ie. not-sky
    pixels).  This flag image can be useful for other things.


    Examples
    --------

    >>> mastersky(data)

    """

    # get the config path
    path = Config().confpath

    # load the WFSS file
    if isinstance(data, str):
        data = WFSS.observed(data)
    elif isinstance(data, WFSS):
        pass
    else:
        raise NotImplementedError
    filename = data.filename

    # check if this file is valid
    if data.instrument == 'ACSWFC':
        skyfile = os.path.join(path, 'instruments', data.config.path,
                               data.config.backgrounds['master'])

    elif data.instrument == 'ACSSBC':
        LOGGER.info('ACS/SBC has no sky background to subtract.')
        return
    elif data.instrument == 'WFC3IR':
        LOGGER.info("WFC3/IR sky background should be done with WFC3_Back_Sub")
        LOGGER.info("by Pirzkal: https://github.com/npirzkal/WFC3_Back_sub")
        return
    elif data.instrument == 'WFC3UVIS':
        LOGGER.info('WFC3/UVIS has no sky background to subtract.')
        return
    elif data.instrument == 'NIRISS':
        raise NotImplementedError(f'{data.instrument}')
    else:
        raise NotImplementedError(f"Instrument {data.instrument} unsupported.")

    # print a message
    LOGGER.info(f"Subtracting master sky for {filename}")

    # how to open the fits files
    if inplace:
        mode = 'update'
    else:
        mode = 'readonly'

    # create some images of the object masks
    phdul = fits.HDUList()
    phdu = fits.PrimaryHDU()
    phdul.append(phdu)

    # process the fits file
    with fits.open(filename, mode=mode) as dhdul,\
            fits.open(skyfile, mode='readonly') as shdul:

        # update the primary header
        dhdul[0].header.set('MSKYSUB', value=True,
                            comment='Flag that sky subtraction done')
        dhdul[0].header.set('MSKYNSG', value=nsig,
                            comment='Number of sigma for source detecting')
        dhdul[0].header.set("MSKYEPS", value=epsilon,
                            comment='Threshold for convergence')
        dhdul[0].header.set("MSKYITER", value=maxiter,
                            comment='Maximum number of iterations')
        dhdul[0].header.set("MSKYFILE", value=skyfile,
                            comment='Name of master sky image')
        headers.add_stanza(dhdul[0].header, 'Master Sky Subtraction',
                           before='MSKYSUB')

        # process each extension
        for hdu in dhdul:

            # grab the name/version of this extension
            name = hdu.header.get("EXTNAME")
            vers = hdu.header.get("EXTVER")

            # only process "SCI" extensions
            if name == "SCI":

                # the science extensions
                sext = ('SCI', vers)

                # read the sky image
                sky = shdul[sext].data
                hdr = shdul[sext].header

                # read the science image
                sci = dhdul[sext].data
                unc = dhdul[("ERR", vers)].data
                dqa = dhdul[("DQ", vers)].data

                # let's be good, and do inv-variance weighting.
                sci2 = sci/unc
                sky2 = sky/unc

                # good pixels based on DQA
                gpx = dqa == 0

                # got to get an initial estimate of the sky level and pixels
                a, m, s = sigma_clipped_stats(sci, sigma=3.,
                                              mask=np.logical_not(gpx))

                # good SKY pixels based on sigma thresholding
                spx = (sci2 - m/unc) < nsig

                # initialize the iteration
                count = 0
                oldcoef = np.inf
                proceed = True

                # start the loop
                while proceed:

                    # find the good pixels
                    g = np.where(np.logical_and(spx, gpx))

                    # extract the good pixels
                    sci3 = sci2[g]
                    sky3 = sky2[g]
                    npix = len(g[0])

                    # do the optimization of a linear function
                    foo = np.sum(sci3*sci3)
                    fom = np.sum(sci3*sky3)
                    fmm = np.sum(sky3*sky3)
                    coef = fom/fmm            # the optimized coefficient
                    chi2 = foo-coef*fom       # the *UNREDUCED* chi2
                    if coef < 0:
                        msg = "Negative sky coefficient encountered.  Exitting."
                        LOGGER.warning(msg)
                        return

                    # update the sky pixel and iteration count
                    spx = (sci2 - coef*sky2) < nsig
                    count += 1

                    # test to stop the loop
                    unconverged = np.abs(oldcoef-coef) > (epsilon*coef)
                    unmaxed = count < maxiter
                    proceed = unconverged and unmaxed

                    # update the state
                    oldcoef = coef

                # save the object mask
                h = hdr.copy()
                h['EXTNAME'] = 'sources'
                phdu = fits.ImageHDU(data=spx.astype(int), header=h)
                phdul.append(phdu)

                # a warning message
                if not unmaxed:
                    LOGGER.warning("Maximum number of iterations reached.")

                res = sci-coef*sky
                src = np.logical_and(spx, gpx)
                sky = np.logical_not(src)

                # print(np.sum(sci*sky)/np.sum(sky))
                # ref=cartesian(sci,sky,filename+f'_{vers}')
                # ldjf

                # update the outputs
                dhdul[sext].data = res  # (sci-coef*sky)

                dhdul[sext].header.set('MSKYSUB', value=True,
                                       comment='Flag that sky subtraction done')
                dhdul[sext].header.set('MSKYCOEF', value=coef,
                                       comment='Scaling coef for master sky')
                dhdul[sext].header.set('MSKYCHI2', value=chi2,
                                       comment='Chi2 of the master sky subtraction')
                dhdul[sext].header.set('MSKYNPIX', value=npix,
                                       comment='Number of pixels used')
                dhdul[sext].header.set('MSKYREDC', value=coef/(npix-1),
                                       comment='Reduced chi2 = chi2/(npix-npar)')
                dhdul[sext].header.set('MSKYITER', value=count,
                                       comment='Number of required iterations')

                headers.add_stanza(dhdul[sext].header, 'Master Sky Subtraction',
                                   before='MSKYSUB')

        # save the object masks
        base = os.path.splitext(os.path.basename(filename))[0]
        phdul.writeto(f'{base}_src.fits', overwrite=True)

        # done with all detectors.  Figure out what to write and  get an
        # output file name
        if inplace:
            outfile = filename
        else:
            if newfile is None:
                # get a new filename
                newfile = f'{base}_msky.fits'

            dhdul.writeto(newfile, overwrite=True)
            outfile = newfile
    return outfile
