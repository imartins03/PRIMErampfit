import numpy as np
from astropy.io import fits
import pandas as pd
import os

# Define directories
input_dir_poly = r'F:\leftover_C1_dif_degrees_test_rampfit\239_frames\unweighted'
input_dir_leg = r'F:\legfit\239_frames_unweighted'
output_dir_residuals = r'F:\legfit\239_frames_unweighted\residuals'
output_csv = r'F:\legfit\239_frames_unweighted\average_residuals_per_frame.csv'

# Number of degrees
num_degrees = 10


def compute_and_save_residuals(degree):
    """
    Compute and save residuals for a given degree.
    """
    # Load coefficient images for the given degree
    file_poly = os.path.join(input_dir_poly, f'fit_coeff_poly_{degree}deg_239frames_noframe1.fits')
    file_leg = os.path.join(input_dir_leg, f'coefficients_leg_239frames_{degree}deg_final_vst_unweighted.fits')

    coeff_poly = fits.getdata(file_poly)
    coeff_leg = fits.getdata(file_leg)

    # Compute the residuals
    residuals = coeff_poly - coeff_leg

    # Define output residual file path
    residual_file = os.path.join(output_dir_residuals, f'residuals_degree_{degree}.fits')

    # Save the residuals to a FITS file
    fits.writeto(residual_file, residuals, overwrite=True)
    print(f"Saved residuals to FITS file: {residual_file}")


def compute_average_residuals():
    """
    Compute the average residuals across all degrees.
    """
    residual_files = [os.path.join(output_dir_residuals, f'residuals_degree_{degree}.fits') for degree in
                      range(1, num_degrees + 1)]

    residuals_list = []
    for residual_file in residual_files:
        with fits.open(residual_file) as hdul:
            residuals = hdul[0].data
        residuals_list.append(residuals)

    # Stack residuals along the degrees dimension (assuming the first dimension is degrees)
    stacked_residuals = np.stack(residuals_list, axis=0)

    # Compute average residuals for each frame
    avg_residuals = np.mean(stacked_residuals, axis=0)  # Averaging over the degrees dimension

    # Compute average residual for each frame
    avg_residual_per_frame = np.mean(avg_residuals, axis=(1, 2))  # Averaging over (1022, 4088) dimensions

    return avg_residual_per_frame


def save_residuals_to_csv(avg_residuals):
    """
    Save average residuals to a CSV file.
    """
    # Create a DataFrame for the CSV file
    frame_numbers = np.arange(1, len(avg_residuals) + 1)  # Frame numbers start from 1
    df = pd.DataFrame({'Frame Number': frame_numbers, 'Average Residual': avg_residuals})

    # Save to CSV
    df.to_csv(output_csv, index=False)
    print(f"Saved residuals to CSV: {output_csv}")


# Create the output directory for residuals if it does not exist
os.makedirs(output_dir_residuals, exist_ok=True)

# Compute and save residuals for degrees from 1 to 10
for degree in range(1, num_degrees + 1):
    compute_and_save_residuals(degree)

# Compute average residuals across all degrees
avg_residuals = compute_average_residuals()

# Save average residuals to CSV
save_residuals_to_csv(avg_residuals)
