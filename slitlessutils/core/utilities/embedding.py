import os
import shutil
from astropy.io import fits
import numpy as np


def embedsub_full_chip(subarray_file, y_axis, x_axis, output_dir=''):
    """
    Embed a subarray flt/flc into a full chip
    The original file is saved as <output_dir>/<root>_original.fits
    Intended for WFC3 and ACS instruments (UVIS, IR, WFC, SBC)

    This function is an adaptation from the G280 transit notebook code:
    https://github.com/spacetelescope/hst_notebooks/blob/main/notebooks/WFC3/uvis_g280_transit/g280_transit_tools.py#L154

    Parameters
    ----------
    subarray_file : str
        Path to original subarray _flt.fits file. e.g. '/my/path/ipppssoot_flt.fits'
    y_axis : int
        Number of y detector rows for a full chip; e.g. 2051 for UVIS, 2048 for WFC, 1014 for IR
    x_axis : int
        Number of x detector rows for a full chip; e.g. 4096 for UVIS and WFC
    output_dir : str
        Optional. Directory to save embedded full-chip file e.g. '/my/other/path/'

    Returns
    -------
    embedded_file : str
        Path to embedded file <output_dir>/<root>.fits
    """

    # Get rootname and build filenames
    filename = os.path.basename(subarray_file)
    root, _ = os.path.splitext(filename)
    dirname = os.path.dirname(subarray_file)

    # Check that subarray is either a flt or flc FITS file
    if not filename.endswith(('_flt.fits', '_flc.fits')):
        raise ValueError("Expected filename ending with '_flt.fits' or '_flc.fits'")

    # Copy the original file to new name
    original_file = os.path.join(output_dir, f"{root}_original.fits")
    shutil.copyfile(subarray_file, original_file)

    # Set full path+filename for the output file
    embedded_file = os.path.join(output_dir, filename)

    # Copy file to output_dir if it's not the cwd
    if output_dir != dirname:
        shutil.copyfile(subarray_file, embedded_file)

    # Open the embedded file and update it
    with fits.open(embedded_file, mode='update') as hdu:

        # Prepare empty full-frame arrays
        sci = np.zeros((y_axis, x_axis), dtype=np.float32)
        err = np.zeros((y_axis, x_axis), dtype=np.float32)
        dq = np.zeros((y_axis, x_axis), dtype=np.int16) + 4

        # Extract header info from SCI extension
        naxis1 = hdu['SCI', 1].header['NAXIS1']
        naxis2 = hdu['SCI', 1].header['NAXIS2']
        ltv1 = hdu['SCI', 1].header['LTV1']
        ltv2 = hdu['SCI', 1].header['LTV2']
        crpix1 = hdu['SCI', 1].header['CRPIX1']
        crpix2 = hdu['SCI', 1].header['CRPIX2']

        x_min = int(-ltv1)
        x_max = x_min + naxis1
        y_min = int(-ltv2)
        y_max = y_min + naxis2

        # Embed the subarray into the full-chip
        sci[y_min:y_max, x_min:x_max] = hdu['SCI', 1].data
        err[y_min:y_max, x_min:x_max] = hdu['ERR', 1].data
        dq[y_min:y_max, x_min:x_max] = hdu['DQ', 1].data

        # Embed samp and time arrays if WFC3/IR subarray
        if hdu[0].header['DETECTOR'] == 'IR':
            samp = np.zeros((x_axis, y_axis), dtype=np.int16)
            time = np.zeros((x_axis, y_axis), dtype=np.float32)
            samp[y_min:y_max, x_min:x_max] = hdu['SAMP', 1].data
            time[y_min:y_max, x_min:x_max] = hdu['TIME', 1].data
            hdu['SAMP', 1].data = samp
            hdu['TIME', 1].data = time

        # Update headers and data
        hdu['SCI', 1].header['SIZAXIS1'] = x_axis
        hdu['SCI', 1].header['SIZAXIS2'] = y_axis

        for ext in ['SCI', 'ERR', 'DQ']:
            if 'CRPIX1' in hdu[ext, 1].header:
                hdu[ext, 1].header['CRPIX1'] = crpix1 + x_min
                hdu[ext, 1].header['CRPIX2'] = crpix2 + y_min
            hdu[ext, 1].header['LTV1'] = 0.0
            hdu[ext, 1].header['LTV2'] = 0.0

        hdu[0].header['SUBARRAY'] = False
        hdu['SCI', 1].data = sci
        hdu['ERR', 1].data = err
        hdu['DQ', 1].data = dq

        hdu[0].header.add_history("This file was modified with function `embedsub_full_chip`")
        hdu[0].header.add_history(f"    LTV1 and LTV2 were originally: {ltv1}, {ltv2}")
        hdu[0].header.add_history(f"    NAXIS1 and NAXIS2 were originally: {naxis1}, {naxis2}")

    return embedded_file


def embedsub_full_detector(subarray_file, y_axis, x_axis, output_dir=''):
    """
    Embeds a full chip subarray into a full detector file, creating blank arrays for unused chip
    Intended for WFC3/UVIS and ACS/WFC

    Parameters
    ----------
    subarray_file : str
        Path to full chip file. e.g. '/my/path/ipppssoot_fit.fits'
    y_axis : int
        Number of y detector rows for a full chip; e.g. 2051 for UVIS, 2048 for WFC, 1014 for IR
    x_axis : int
        Number of x detector rows for a full chip; e.g. 4096 for UVIS & WFC
    output_dir : str
        Optional. Directory to save embedded full-detector file e.g. '/my/path/dir'

    Returns
    -------
    embedded_file : str
        Path to the embedded full-frame file <output_dir>/<root>.fits
    """

    # Embed subarray into full chip
    embedded_file = embedsub_full_chip(subarray_file, y_axis, x_axis, output_dir)

    # Create full detector FITS
    with fits.open(embedded_file, mode='update') as hdu:

        # Grab the real SCI/ERR/DQ (EXTVER=1)
        sci1 = hdu['SCI', 1]
        err1 = hdu['ERR', 1]
        dq1 = hdu['DQ', 1]

        # Determine the blank chip
        ccdchip_real = sci1.header.get('CCDCHIP')
        ccdchip_blank = 2 if ccdchip_real == 1 else 1

        # Create blank HDUs
        sci_blank = sci1.copy()
        err_blank = err1.copy()
        dq_blank = dq1.copy()

        sci_blank.header['CCDCHIP'] = ccdchip_blank
        sci_blank.data = np.zeros_like(sci1.data)
        err_blank.data = np.zeros_like(err1.data)
        dq_blank.data = np.zeros_like(dq1.data, dtype=np.int16) + 4

        # Assign EXTVERs according to real chip
        if ccdchip_real == 1:
            # Blank CCDCHIP=2 -> EXTVER=1
            sci_blank.header['EXTVER'] = 1
            err_blank.header['EXTVER'] = 1
            dq_blank.header['EXTVER'] = 1

            # Real CCDCHIP=1 -> EXTVER=2
            sci1.header['EXTVER'] = 2
            err1.header['EXTVER'] = 2
            dq1.header['EXTVER'] = 2

            # Insert blank HDUs before real HDUs
            hdu.insert(1, dq_blank)
            hdu.insert(1, err_blank)
            hdu.insert(1, sci_blank)

        else:
            # Real CCDCHIP=2 -> EXTVER=1
            sci1.header['EXTVER'] = 1
            err1.header['EXTVER'] = 1
            dq1.header['EXTVER'] = 1

            # Blank CCDCHIP=1 -> EXTVER=2
            sci_blank.header['EXTVER'] = 2
            err_blank.header['EXTVER'] = 2
            dq_blank.header['EXTVER'] = 2

            # Insert blank HDUs after real HDUs
            hdu.insert(4, dq_blank)
            hdu.insert(4, err_blank)
            hdu.insert(4, sci_blank)

        # Save changes in place
        hdu.flush()

    return embedded_file
