from astropy.io import fits
import numpy as np
from astropy.stats import sigma_clipped_stats
from photutils.detection import find_peaks
from photutils.segmentation import detect_sources
from scipy.ndimage import median_filter

def reduce_science_frame(
    science_filename,
    median_bias_filename,
    median_flat_filename,
    median_dark_filename,
    reduced_science_filename="reduced_science.fits",
):
    """Reduce a science frame using bias, dark, and flat frames."""

    # load all data
    with fits.open(science_filename) as hdul:
        science_header = hdul[0].header
        exposure_time = science_header.get('EXPTIME')
        science_data = fits.getdata(science_filename).astype('f4')

    bias_data = fits.getdata(median_bias_filename).astype('f4')
    dark_data = fits.getdata(median_dark_filename).astype('f4')
    flat_data = fits.getdata(median_flat_filename).astype('f4')
    

    # subtract bias
    reduced_science = science_data - bias_data

    # subtract dark (scaled by exposure time)
    reduced_science -= dark_data * exposure_time

    # flat field correction
    reduced_science /= flat_data

    # remove cosmic rays (simple median filter + threshold method)
    # this can be replaced with astroscrappy or similar for better accuracy
    smoothed = median_filter(reduced_science, size=3)
    residual = reduced_science - smoothed
    threshold = 5 * np.std(residual)
    mask = residual > threshold
    reduced_science[mask] = smoothed[mask]  # Replace cosmic rays

    # save to FITS
    primary_science = fits.PrimaryHDU(data=reduced_science.data, header=fits.Header())
    primary_science.writeto(reduced_science_filename, overwrite=True)

    return reduced_science
