from astropy.io import fits
import numpy as np

def calculate_gain(files):
    """
    Calculate the detector gain (e-/ADU) using two flat-field images.

    Parameters:
        files (list): List of two flat-field FITS file paths.

    Returns:
        float: Calculated gain in e-/ADU.
    """
    if len(files) != 2:
        raise ValueError("You must provide exactly two flat-field images.")

    # Read the flat-field images
    flat1 = fits.getdata(files[0]).astype('f4')
    flat2 = fits.getdata(files[1]).astype('f4')

    # Calculate mean and variance
    mean1 = np.mean(flat1)
    mean2 = np.mean(flat2)
    mean_comb = (mean1 + mean2) / 2.0

    diff = flat1 - flat2
    variance_diff = np.var(diff) / 2.0  # Because var(A - B) = 2 * var(single image)

    gain = mean_comb / variance_diff

    return float(gain)

def calculate_readout_noise(files, gain):
    """
    Calculate the readout noise in electrons (e-) using two bias frames.

    Parameters:
        files (list): List of two bias FITS file paths.
        gain (float): Gain in e-/ADU.

    Returns:
        float: Readout noise in electrons.
    """
    if len(files) != 2:
        raise ValueError("You must provide exactly two bias frames.")

    # Read the bias frames
    bias1 = fits.getdata(files[0]).astype(np.float32)
    bias2 = fits.getdata(files[1]).astype(np.float32)

    # Calculate variance of the difference
    diff = bias1 - bias2
    variance_diff = np.var(diff) / 2.0  # Same reasoning as above

    readout_noise = np.sqrt(variance_diff) * gain

    return float(readout_noise)
