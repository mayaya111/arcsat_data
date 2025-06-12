import numpy as np
from astropy.io import fits 
from astropy.stats import sigma_clip
import pathlib
import os

def create_median_bias(bias_list, median_bias_filename):
    """Creates a median bias frame from a list of FITS bias files.

    Parameters:
    - bias_list: list of strings, paths to bias FITS files.
    - median_bias_filename: string, path to save the resulting median bias frame.

    Returns:
    - A 2D numpy array of the median bias frame.
    """
    # read all bias frames into a 3D numpy array (stack of 2D images)
    bias_frames = []
    for file in bias_list:
        data = fits.getdata(file).astype('f4')
        bias_frames.append(data)
        
    bias_stack = np.array(bias_frames)
    
    # apply sigma clipping along the stack axis (axis=0)
    clipping = sigma_clip(bias_stack, cenfunc = 'median', sigma=3.0, axis=0)
    

    # take the median of the clipped data
    # convert masked array to regular ndarray before saving
    median_bias = np.ma.median(clipping, axis=0)
    if isinstance(median_bias, np.ma.MaskedArray):
        median_bias = median_bias.filled()
    # returns fits immages 
    primary = fits.PrimaryHDU(data=median_bias, header=fits.Header())
    hdul = fits.HDUList([primary])
    hdul.writeto(median_bias_filename, overwrite=True)

    return median_bias