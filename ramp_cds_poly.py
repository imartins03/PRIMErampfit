# import numpy as np
# from astropy.io import fits
# import os
#
# # File paths
# first_filename = r'D:\NLC\C1\01124973C1.fits.fz'  # First after superbias
# last_num = 1124973 + 293
# last_filename = f'D:\\NLC\\C1\\0{last_num}C1.fits.fz'  # Corrected string formatting
#
# n_frames = last_num - 1124973
#
# # Load FITS data
# first = fits.getdata(first_filename)
# last = fits.getdata(last_filename)
#
#
# def get_ramp_cds():
#     # Calculate CDS (Correlated Double Sampling)
#     cds = last - first
#
#     # Define output filename
#     output_filename = f'D:\\NLC\\C1\\last-first_{n_frames}frames.cds.fits'
#
#     # Write to FITS file
#     fits.writeto(output_filename, cds, overwrite=True)
#     print(f'CDS file written to: {output_filename}')
#
#
# get_ramp_cds()
# print(last_filename)

import numpy as np
from astropy.io import fits
import os

# # Define file path template
# fit_cube_path_template = r'D:\NLC\C1\dif_degrees_test\100_frames\fit_cube_poly_{degree}deg_100frames_noframe1.fits'
# degree = 6
#
# # Construct the FITS file path
# fit_cube_path = fit_cube_path_template.format(degree=degree)
#
# def subtract_first_last_frame(fit_cube_path):
#     # Load FITS cube data
#     fit_cube = fits.getdata(fit_cube_path)
#
#     # Extract the first and last frames
#     first_frame = fit_cube[0]
#     last_frame = fit_cube[-1]
#
#     # Calculate the difference (last frame - first frame)
#     difference = last_frame - first_frame
#
#     # Define output filename
#     output_filename = fit_cube_path.replace('fit_cube_poly', 'last_minus_first_frame')
#
#     # Write the result to a new FITS file
#     fits.writeto(output_filename, difference, overwrite=True)
#     print(f'Difference file written to: {output_filename}')
#
# # Call the function with the FITS cube path
# subtract_first_last_frame(fit_cube_path)

file_for_239_raw = r"D:\NLC\C1\last-first_293frames.cds.fits"
file_for_100_raw = r"D:\NLC\C1\last-first_100frames.cds.fits"

file_for_239_fit = r"F:\leftover_C1_dif_degrees_test_rampfit\239_frames\last_minus_first_frame_6deg_239frames_noframe1.fits"
file_for_100_fit = r"D:\NLC\C1\dif_degrees_test\100_frames\last_minus_first_frame_6deg_100frames_noframe1.fits"

file_for_239_raw = r"D:\NLC\C1\last-first_293frames.cds.fits"
file_for_100_raw = r"D:\NLC\C1\last-first_100frames.cds.fits"

# File paths for fit frames
file_for_239_fit = r"F:\leftover_C1_dif_degrees_test_rampfit\239_frames\last_minus_first_frame_6deg_239frames_noframe1.fits"
file_for_100_fit = r"D:\NLC\C1\dif_degrees_test\100_frames\last_minus_first_frame_6deg_100frames_noframe1.fits"


def compute_ratios():
    # Load raw CDS data
    data_239_raw = fits.getdata(file_for_239_raw)
    data_100_raw = fits.getdata(file_for_100_raw)

    # Calculate the ratio for raw CDS data
    ratio_raw = np.divide(data_239_raw, data_100_raw)

    # Define output filename for raw ratio
    ratio_file_raw = r"D:\NLC\C1\ratio_239_to_100_raw.fits"

    # Save the ratio for raw data
    fits.writeto(ratio_file_raw, ratio_raw, overwrite=True)
    print(f'Raw CDS ratio file written to: {ratio_file_raw}')

    # Load fit frame data
    data_239_fit = fits.getdata(file_for_239_fit)
    data_100_fit = fits.getdata(file_for_100_fit)

    # Calculate the ratio for fit data
    ratio_fit = np.divide(data_239_fit, data_100_fit)

    # Define output filename for fit ratio
    ratio_file_fit = r"D:\NLC\C1\ratio_239_to_100_fit.fits"

    # Save the ratio for fit data
    fits.writeto(ratio_file_fit, ratio_fit, overwrite=True)
    print(f'Fit ratio file written to: {ratio_file_fit}')

    return ratio_file_raw, ratio_file_fit


# Call the function to compute ratios and save files
compute_ratios()
