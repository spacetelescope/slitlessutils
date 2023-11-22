import os

import numpy as np
from astropy.io import fits
from scipy import ndimage
from skimage import measure, morphology

from ....logger import LOGGER
from ...utilities import headers, indices
from .crbitvalues import BITVALUES

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


def _laplace(filename, inplace=True, newfile=None, bitvalue=None,
             kernel='3a', snlim=9., minpix=3, maxpix=100,
             growsize=2, growform='diamond'):
    """
    Helper function to identify cosmic rays (CRs) by Laplace filtering.

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
        # if bitvalue is an int, check if its a multiple of 2
        if isinstance(bitvalue, int):
            get_header_bitvalue = not np.log2(bitvalue).is_integer()
        else:
            get_header_bitvalue = True

        # if bitvalue is invalid, get something from the header
        if get_header_bitvalue:

            # grab the instrument name
            ins = hdul[0].header.get('INSTRUME')

            # grab the bitvalue to mask
            bitvalue = BITVALUES[ins]

        # process each extension
        for hdu in hdul:
            if hdu.header.get('EXTNAME', 'primary') == 'SCI':

                # get the extension version for use below
                extver = hdu.header.get('EXTVER', 1)

                # convolve with a Laplacian and weight by uncertainty
                sci = hdul[('SCI', extver)].data
                con = ndimage.convolve(sci, kern)
                s2n = np.abs(con) / hdul[('ERR', extver)].data

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
            # hdul[0].header.set('CRINTERP', value=interp,
            #                   comment='Interpolate over CRs in science exten')
            # hdul[0].header.set('CRINTMTH', value=interpmethod,
            #                   comment='method for CR interpolation')
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


def laplace(filenames, **kwargs):
    """
    Function to perform the Laplace filtering of cosmic rays

    Parameters
    ----------
    filenames : str or list/tuple
        The name of files to search for cosmic rays

    **kwargs : dict
        See `_laplace()`

    Return
    ------
    outfiles : str or list/tuple
        The name of the files with the cosmic rays masked
    """

    if isinstance(filenames, (list, tuple)):
        outfiles = [_laplace(f, **kwargs) for f in filenames]
    elif isinstance(filenames, str):
        outfiles = _laplace(filenames, **kwargs)
    else:
        LOGGER.warning(f'Invalid data type for filenames: {type(filenames)}')

    return outfiles


# if __name__=='__main__':
#
#    filename='/Users/rryan/ACS/WFC/HST/JDN603020/jdn603drq_flc.fits'
#    filename='jdn603drq_flc.fits'
#    x=laplace(filename)
