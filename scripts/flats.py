from astropy.io import fits
import numpy as np
import os
from collections import defaultdict

def create_median_flats_per_filter(flat_files, bias_file, dark_file=None, output_dir="."):
    """
    Create median flats grouped by filter.

    flat_files: list of flat FITS file paths
    bias_file: path to median bias FITS
    dark_file: optional median dark FITS path
    output_dir: where to save median flats
    
    Saves median flats named median_flat_<filter>.fits
    """
    # load bias and dark data once
    bias_data = fits.getdata(bias_file)
    dark_data = fits.getdata(dark_file) if dark_file else None
    
    # group flats by FILTER header
    flats_by_filter = defaultdict(list)
    for f in flat_files:
        with fits.open(f) as hdul:
            filter_name = hdul[0].header.get("FILTER", "UNKNOWN")
        flats_by_filter[filter_name].append(f)
    
    for filter_name, files in flats_by_filter.items():
        print(f"Processing {len(files)} flats for filter {filter_name}")
        calibrated_flats = []
        for fname in files:
            flat_data = fits.getdata(fname)
            # Bias subtract
            flat_corr = flat_data - bias_data
            # Dark subtract if available
            if dark_data is not None:
                flat_corr -= dark_data
            calibrated_flats.append(flat_corr)
        
        # stack and median combine
        stack = np.array(calibrated_flats)
        median_flat = np.median(stack, axis=0)
        
        # normalize
        median_val = np.median(median_flat)
        if median_val == 0:
            raise ValueError(f"Median of flat for filter {filter_name} is zero")
        median_flat /= median_val
        
        # save
        outname = os.path.join(output_dir, f"median_flat_{filter_name}.fits")
        fits.PrimaryHDU(data=median_flat).writeto(outname, overwrite=True)
        print(f"Saved median flat for filter {filter_name} to {outname}")

