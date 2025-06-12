#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @Filename: darks.py
# @License: BSD 3-clause (http://www.opensource.org/licenses/BSD-3-Clause)


from astropy.io import fits
from astropy.stats import sigma_clip
import numpy as np

def create_median_dark(dark_list, bias_filename, median_dark_filename):
    """Creates a median dark frame corrected for bias and normalized by exposure time.

    Parameters:
    - dark_list: list of strings, paths to dark FITS files.
    - bias_filename: string, path to the median bias FITS file.
    - median_dark_filename: string, path to save the resulting median dark frame.

    Returns:
    - A 2D numpy array of the median dark frame.
    """

    # load the median bias frame
    bias_data = fits.getdata(bias_filename).astype('f4')

    # list to store bias-subtracted and normalized dark frames
    dark_corrected_frames = []

    for path in dark_list:
        with fits.open(path) as dark:
            dark_data = dark[0].data.astype('f4')
            exptime = dark[0].header['EXPTIME']
            if exptime <= 0:
                raise ValueError(f"Invalid EXPTIME ({exptime}) in {path}")
            dark_nobias = dark_data - bias_data
            dark_corrected_frames.append(dark_nobias / exptime)
            header = dark[0].header.copy()


    # stack and sigma clip the corrected dark frames
    dark_3d = np.array(dark_corrected_frames)
    clipping = sigma_clip(dark_3d, cenfunc = 'median', sigma=3.0, axis=0)

    # compute the median of the clipped data
    median_dark = np.ma.median(clipping, axis=0).filled()

    # save to a FITS file
    dark_hdu = fits.PrimaryHDU(data=median_dark, header=fits.Header())
    dark_hdu.writeto(median_dark_filename, overwrite=True)

    return median_dark
