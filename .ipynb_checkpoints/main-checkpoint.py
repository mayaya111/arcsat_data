import os
from pathlib import Path
from astropy.io import fits
import numpy as np
from photutils.aperture import CircularAperture, CircularAnnulus, aperture_photometry

from ccd.bias import create_median_bias
from ccd.darks import create_median_dark
from ccd.flats import create_median_flat
from ccd.science import reduce_science_frame


def run_reduction(data_dir):
    data_dir = Path(data_dir)
    print("Starting CCD Reduction Pipeline...\n")

    # Find files
    bias_list = sorted([str(f) for f in data_dir.glob('Bias-*.fit')])
    dark_list = sorted([str(f) for f in data_dir.glob('Dark-*.fit')])
    flat_list = sorted([str(f) for f in data_dir.glob('AutoFlat-*.fit')])
    science_files = sorted([str(f) for f in data_dir.glob('kelt-16-*.fit')])

    print(f"Found {len(bias_list)} bias files")
    print(f"Found {len(dark_list)} dark files")
    print(f"Found {len(flat_list)} flat files")
    print(f"Found {len(science_files)} science files")

    # Create calibration frames
    print("Creating median bias...")
    median_bias_path = data_dir / "median_bias.fits"
    create_median_bias(bias_list, str(median_bias_path))

    print("Creating median dark...")
    median_dark_path = data_dir / "median_dark.fits"
    create_median_dark(dark_list, str(median_bias_path), str(median_dark_path))

    print("Creating median flat...")
    median_flat_path = data_dir / "median_flat.fits"
    create_median_flat(flat_list, str(median_bias_path), str(median_flat_path), dark_filename=str(median_dark_path))

    # Reduce science frames
    reduced_files = []
    for sci_file in science_files:
        outname = data_dir / ("reduced_" + os.path.basename(sci_file))
        reduce_science_frame(
            science_filename=sci_file,
            median_bias_filename=str(median_bias_path),
            median_dark_filename=str(median_dark_path),
            median_flat_filename=str(median_flat_path),
            reduced_science_filename=str(outname)
        )
        print(f"Reduced science frame saved to {outname}")
        reduced_files.append(str(outname))

    # Optional: Perform aperture photometry on first reduced science frame
    if reduced_files:
        first_file = reduced_files[0]
        with fits.open(first_file) as hdul:
            data = hdul[0].data

        positions = [(100, 100)]  # example position, adjust as needed
        aperture_radius = 5
        annulus_r_in = 8
        annulus_r_out = 12

        apertures = CircularAperture(positions, r=aperture_radius)
        annuli = CircularAnnulus(positions, r_in=annulus_r_in, r_out=annulus_r_out)

        bkg_medians = []
        for annulus in annuli:
            mask = annulus.to_mask(method='exact')
            annulus_data = mask.multiply(data)
            annulus_data_1d = annulus_data[mask.data > 0]
            bkg_medians.append(np.median(annulus_data_1d))

        phot_table = aperture_photometry(data, apertures)
        for i in range(len(phot_table)):
            phot_table['aperture_sum'][i] -= bkg_medians[i] * apertures.area

        print("\nAperture Photometry Results:")
        print(phot_table)

    print("\nCCD Reduction Pipeline complete.")
