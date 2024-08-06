import numpy as np
from astropy.io import fits
import os

# File paths
first_filename = r'D:\NLC\C1\01124973C1.fits.fz'  # First after superbias
last_num = 1124973 + 293
last_filename = f'D:\\NLC\\C1\\0{last_num}C1.fits.fz'  # Corrected string formatting

n_frames = last_num - 1124973

# Load FITS data
first = fits.getdata(first_filename)
last = fits.getdata(last_filename)


def get_ramp_cds():
    # Calculate CDS (Correlated Double Sampling)
    cds = last - first

    # Define output filename
    output_filename = f'D:\\NLC\\C1\\last-first_{n_frames}frames.cds.fits'

    # Write to FITS file
    fits.writeto(output_filename, cds, overwrite=True)
    print(f'CDS file written to: {output_filename}')


get_ramp_cds()
print(last_filename)