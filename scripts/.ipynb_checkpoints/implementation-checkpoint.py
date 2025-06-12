
from pathlib import Path
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# define the data directory
data_dir = Path("astr-480-env/work/arcsat_data/scripts/data")  # Adjust path depending on location of .fit files

# load and ort your calibration frames
bias_files = sorted([str(f) for f in data_dir.glob("Bias_*.fits")])
dark_files = sorted([str(f) for f in data_dir.glob("Dark_*.fits")])
flat_files = sorted([str(f) for f in data_dir.glob("domeflat_*.fits")])
sci_files = sorted([str(f) for f in data_dir.glob("Cats_*.fits")])
print("Bias files:", bias_files)

#implementing functions
create_median_bias(bias_files, str(data_dir / "median_bias.fits"))
create_median_dark(dark_files, str(data_dir / "median_bias.fits"), str(data_dir / "median_dark.fits"))
create_median_flats_per_filter(
    flat_files=flat_files,
    bias_file=str(median_bias_file),
    dark_file=str(median_dark_file),
    output_dir=str(data_dir)
)

#making sure the median flats and images are sorted in the same way 
median_flats = sorted(data_dir.glob("median_flat_*.fits"))
images = sorted(data_dir.glob("Cats_*.fits"))

# reduce science frame manually for each science frame 
reduce_science_frame(
    science_filename=sci_files[0],
    median_bias_filename=str(data_dir / "median_bias.fits"),
    median_dark_filename=str(data_dir / "median_dark.fits"),
    median_flat_filename=str(median_flats[0]),
    reduced_science_filename = str(data_dir / ("reduced_" + Path(sci_files[0]).name))
)

reduce_science_frame(
    science_filename=sci_files[1],
    median_bias_filename=str(data_dir / "median_bias.fits"),
    median_dark_filename=str(data_dir / "median_dark.fits"),
    median_flat_filename=str(median_flats[1]),
    reduced_science_filename = str(data_dir / ("reduced_" + Path(sci_files[1]).name))
)

reduce_science_frame(
    science_filename=sci_files[2],
    median_bias_filename=str(data_dir / "median_bias.fits"),
    median_dark_filename=str(data_dir / "median_dark.fits"),
    median_flat_filename=str(median_flats[2]),
    reduced_science_filename = str(data_dir / ("reduced_" + Path(sci_files[2]).name))
)

reduce_science_frame(
    science_filename=sci_files[3],
    median_bias_filename=str(data_dir / "median_bias.fits"),
    median_dark_filename=str(data_dir / "median_dark.fits"),
    median_flat_filename=str(median_flats[3]),
    reduced_science_filename = str(data_dir / ("reduced_" + Path(sci_files[3]).name))
)

reduce_science_frame(
    science_filename=sci_files[4],
    median_bias_filename=str(data_dir / "median_bias.fits"),
    median_dark_filename=str(data_dir / "median_dark.fits"),
    median_flat_filename=str(median_flats[4]),
    reduced_science_filename = str(data_dir / ("reduced_" + Path(sci_files[4]).name))
)


# show results (wanted to make sure it was running properly so I ran for H-Alpha and r filters
with fits.open(data_dir / f"reduced_{Path(sci_files[0]).name}") as hdul:
    plt.imshow(hdul[0].data, cmap='gray', origin='lower', vmin=0, vmax=np.percentile(hdul[0].data, 99))
    plt.title("Reduced Science Frame")
    plt.colorbar()
    plt.show()

with fits.open(data_dir / f"reduced_{Path(sci_files[3]).name}") as hdul:
    plt.imshow(hdul[0].data, cmap='gray', origin='lower', vmin=0, vmax=np.percentile(hdul[0].data, 99))
    plt.title("Reduced Science Frame")
    plt.colorbar()
    plt.show()
