import numpy as np
from astropy.io import fits
import pandas as pd
import os

# Define directories
input_dir_poly = r'F:\leftover_C1_dif_degrees_test_rampfit\239_frames\unweighted'
input_dir_leg = r'F:\legfit\239_frames_unweighted'
output_dir_residuals = r'F:\legfit\239_frames_unweighted\residuals'
output_csv = r'F:\legfit\239_frames_unweighted\average_residuals_per_frame.csv'

# Degree to process
degree = 1


def compute_and_save_residuals_for_first_degree():
    """
    Compute and save residuals for the first degree.
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


def compute_average_residuals_per_frame():
    """
    Compute the average residuals for each frame across all frames.
    """
    residual_file = os.path.join(output_dir_residuals, f'residuals_degree_{degree}.fits')

    # Load the residuals
    with fits.open(residual_file) as hdul:
        residuals = hdul[0].data

    # Compute the average residual for each frame
    avg_residuals = np.mean(residuals, axis=(1, 2))  # Averaging over (1022, 4088) dimensions

    return avg_residuals


def save_residuals_to_csv(avg_residuals):
    """
    Save average residuals per frame to a CSV file.
    """
    # Create a DataFrame for the CSV file
    frame_numbers = np.arange(1, len(avg_residuals) + 1)  # Frame numbers start from 1
    df = pd.DataFrame({'Frame Number': frame_numbers, 'Average Residual': avg_residuals})

    # Save to CSV
    df.to_csv(output_csv, index=False)
    print(f"Saved average residuals to CSV: {output_csv}")


# Create the output directory for residuals if it does not exist
os.makedirs(output_dir_residuals, exist_ok=True)

# Compute and save residuals for the first degree
compute_and_save_residuals_for_first_degree()

# Compute average residuals for each frame
avg_residuals = compute_average_residuals_per_frame()

# Save average residuals to CSV
save_residuals_to_csv(avg_residuals)
