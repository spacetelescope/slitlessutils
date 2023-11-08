import os

import numpy as np
from astropy.io import fits
from scipy import interpolate, ndimage
from skimage import measure, morphology

from ....logger import LOGGER
from ...utilities import headers, indices

# some pre-built Laplace kernels
KERNELS = {'3a': [[0, -1, 0],
                  [-1, 4, -1],
                  [0, -1, 0]],
           '3b': [[-1, -1, -1],
                  [-1, 8, -1],
                  [-1, -1, -1]],
           '3c': [[1, -2, 1],
                  [-2, 4, -2],
                  [1, -2, 1]],
           '5a': [[0, 0, -1, 0, 0],
                  [0, -1, -2, -1, 0],
                  [-1, -2, 16, -2, -1],
                  [0, -1, -2, -1, 0],
                  [0, 0, -1, 0, 0]]}

# what bit value to put in the DQAs
BITVALUES = {'ACS': 16384}
DEFAULT = 8192


def laplace(filename, inplace=True, newfile=None, bitvalue=None,
            interp=False, interpmethod='cubic',
            kernel='3a', snlim=9., minpix=3, maxpix=100,
            growsize=2, growform='diamond'):
    """
    Function to identify cosmic rays (CRs) by Laplace filtering.

    The default nature is to update the data-quality array (DQA) of the
    file in place, without altering the science or error arrays.

    Parameters
    ----------
    filename : str
        The full path to an 'flt'-like file.

    inplace : bool, optional
        Flag that the updates to the DQA should be done inplace and not
        write a new file.  Default is True

    newfile : str or None, optional
        If the inplace flag is set to False, then this will be the name
        of the updated file.  If newfile is set to None, then the new
        file will have a name of:

        f'{os.path.basename(filename)}_crj.fits'

    bitvalue : int or None, optional
        The bit value to add to the DQA to indicate a CR was identified
        If set to None, then the primary header is searched for a
           keyword 'INSTRUME' which is used with a dictionary of defaults
           hardcoded for different instruments.
        If set to a valid integer (ie. one that is a multiple of 2), then
           that is used.

    interp : bool, optional
        Flag to interpolate over CRs in the science image.  This is note
        advised for anything other than visualization purposes as the
        bad pixels will be rejected and this is a slow operation.  Therefore
        there is little value, other than display. Default is False

    interpmethod : str, optional
        Keyword governing the form of interpolation sent to `np.griddata()`
        Default is 'linear'

    kernel : str, optional
        The form of the Laplacian kernel.  The higher order kernels provide
        a better approximation of the Laplace operator, but come at
        CPU performance and possible increased noise.  Valid options
        are '3a', '3b', '3c', '5a'.  The number represents the size of the
        kernel.  Default is '3a'.

    snlim : float, optional
        The limit on the S/N for a CR.  Default is 9.

    minpix : int, optional
        The minimum size of a CR (in pixels)  Default is 3

    maxpix : int, optional
        The maximum size of a CR (in pixels)  Default is 100

    growize : int or tuple (depending on `growform`), optional
        Size to grow the CR mask (in pixels). For some growing elements
        this is a scalar integer, and for others it is a 2-tuple as shown:

        Grow Form  |  Grow Size
        -----------+-----------
         square    |  scalar
         rectangle |  tuple
         diamond   |  scalar
         disk      |  scalar
         octagon   |  tuple
         star      |  scalar

        Default is 2.

    growform : str, optional
        The form of the growing array.  Valid options are 'square',
        'rectangle', 'diamond', 'disk', 'octagon', and 'star'.
        Default is 'diamond'

    Returns
    -------
    filename : str
        The name of the updated file.

    """

    LOGGER.info(f'Flagging cosmic rays with Laplacian kernel {filename}')

    # growing structure
    grower = None
    if growsize is not None and growsize > 0 and growform is not None:
        growform = growform.lower()
        if growform == 'square':
            grower = morphology.square(growsize)
        elif growform == 'rectangle':
            grower = morphology.rectangle(*growsize)
        elif growform == 'diamond':
            grower = morphology.diamond(growsize)
        elif growform == 'disk':
            grower = morphology.disk(growsize)
        elif growform == 'octagon':
            grower = morphology.octagon(*growsize)
        elif growform == 'star':
            grower = morphology.star(growsize)
        else:
            LOGGER.warning('No valid grow form found.')

    # grab the kernels
    kern = np.array(KERNELS.get(kernel, '3a'), dtype=float)
    if np.sum(kern) != 0:
        LOGGER.error("Invalid kernel: sums to non-zero")
        return

    # normalize the kernel to keep S/N meaningful
    kern /= np.amax(kern)

    # set some flags
    didsomething = False
    if inplace:
        mode = 'update'
    else:
        mode = 'readonly'

    # open the file
    with fits.open(filename, mode=mode) as hdul:

        # assume we have to look in the header to get the INSTRUME to
        # lookup the bitvalue.  This could change if a bitvalue is passed
        # and it is of valid type
        get_header_bitvalue = True

        # if bitvalue is an int, check if its a multple of 2
        if isinstance(bitvalue, int):
            get_header_bitvalue = not np.log2(bitvalue).is_integer()

        # if bitvalue is invalid, get something from the header
        if get_header_bitvalue:

            # grab the instrument name
            ins = hdul[0].header.get('INSTRUME')

            # could do this in a single line: ins = BITVALUES.get(ins,DEFAULT)
            # but I want to log a warning message if using the default, so
            # break it up into an if/else block.
            if ins in BITVALUES:
                # if ins is in the default dict, then use it
                bitvalue = BITVALUES[ins]
            else:
                # else use the actual default
                LOGGER.warning(f"Instrument {ins} not found, using default bit value: {DEFAULT}")
                bitvalue = DEFAULT

        # process each extension
        for hdu in hdul:
            if hdu.header.get('EXTNAME', 'primary') == 'SCI':

                # get the extension version for use below
                extver = hdu.header['EXTVER']

                # convolve with a Laplacian and weight by uncertainty
                sci = hdul[('SCI', extver)].data
                con = ndimage.convolve(sci, kern)
                s2n = np.abs(con)/hdul[('ERR', extver)].data

                # consider zooming?
                # zoom=3.
                # dat=ndimage.zoom(hdul[('SCI',extver)].data,zoom,order=1)
                # err=ndimage.zoom(hdul[('ERR',extver)].data,zoom,order=1)
                # con=ndimage.convolve(dat,kern)
                # s2n=np.abs(con)/err
                # s2n=ndimage.zoom(s2n,1./zoom,order=1)

                # flag based on signal-to-noise limit
                bpx = (s2n > snlim)

                # close for small holes in cosmic ray masks
                # bpx=morphology.binary_closing(bpx,footprint=close)

                # update the bpx with small/large sources
                minpix = 0 if minpix is None else minpix
                maxpix = np.inf if maxpix is None else maxpix
                crs = indices.reverse(measure.label(bpx))
                for i, pix in crs.items():
                    if i > 0:
                        n = len(pix[0])
                        if n < minpix or n > maxpix:
                            bpx[pix] = False

                # grow the final mask
                if grower is not None:
                    bpx = morphology.dilation(bpx, footprint=grower)

                # update the data
                g = np.where(bpx)
                if g[0].size > 0:
                    # interpolate over in the science image?
                    # this should be consider cosmetic *ONLY*
                    if interp:
                        LOGGER.info(f"Interpolating over CRs: {filename}[('SCI',{extver})]")
                        gpx = ~bpx

                        xx, yy = np.meshgrid(np.arange(sci.shape[1], dtype=int),
                                             np.arange(sci.shape[0], dtype=int))

                        crs = indices.reverse(measure.label(bpx))
                        for idx, (y, x) in crs.items():
                            if idx > 0:
                                y0 = max(np.amin(y)-3, 0)
                                x0 = max(np.amin(x)-3, 0)
                                y1 = min(np.amax(y)+3, sci.shape[0]-1)
                                x1 = min(np.amax(x)+3, sci.shape[1]-1)

                                subsci = sci[y0:y1, x0:x1]
                                subgpx = gpx[y0:y1, x0:x1]
                                subbpx = bpx[y0:y1, x0:x1]
                                subx = xx[y0:y1, x0:x1]
                                suby = yy[y0:y1, x0:x1]
                                by, bx = np.where(subbpx)

                                spline = interpolate.SmoothBivariateSpline(subx[subgpx],
                                                                           suby[subgpx],
                                                                           subsci[subgpx],
                                                                           s=1e6, kx=1, ky=1)

                                hdul[('SCI', extver)].data[by+y0, bx+x0] = spline(subx[subbpx],
                                                                                  suby[subbpx],
                                                                                  grid=False)

                        # gxy=(xx[gpx],yy[gpx])
                        # gsci=sci[gpx]
                        # bxy=(xx[bpx],yy[bpx])
                        #
                        # newval=interpolate.griddata(gxy,gsci,bxy,
                        #                            fill_value=np.nan,
                        #                            method=interpmethod)
                        # hdul[('SCI',extver)].data[g]=newval

                    # update the DQ array
                    hdul[('DQ', extver)].data[g] += bitvalue
                    didsomething = True      # update flag that things happened

        if didsomething:

            # update the primary header
            hdul[0].header['HISTORY'] = 'updated DQA for CRs by "laplace.py"'
            hdul[0].header.set('CRKERN', value=kernel, comment='Laplace kernel')
            hdul[0].header.set('CRSN', value=snlim, comment='S/N per pix')
            hdul[0].header.set('CRMINPIX', value=minpix,
                               comment='Minimum CR size')
            hdul[0].header.set('CRINTERP', value=interp,
                               comment='Interpolate over CRs in science exten')
            hdul[0].header.set('CRINTMTH', value=interpmethod,
                               comment='method for CR interpolation')
            hdul[0].header.set('CRGRWSZ', value=str(growsize),
                               comment='Npix to grow CRs')
            hdul[0].header.set('CRGRWFRM', value=growform,
                               comment='Shape to grow')
            headers.add_stanza(hdul[0].header, 'Laplace CR parameters',
                               before='CRKERN')

        # juggle the outputting
        if inplace:
            outfile = filename

        else:
            if newfile is None:
                # if a newfilename is not given parse the input name
                base = os.path.splitext(os.path.basename(filename))[0]
                newfile = f'{base}_crj.fits'
            hdul.writeto(newfile, overwrite=True)
            outfile = newfile

    return outfile


# if __name__=='__main__':
#
#    filename='/Users/rryan/ACS/WFC/HST/JDN603020/jdn603drq_flc.fits'
#    filename='jdn603drq_flc.fits'
#    x=laplace(filename)
