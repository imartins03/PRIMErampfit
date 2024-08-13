from astropy.io import fits
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Define paths
y_cube_path = r'D:\NLC\C1\y_cube_500.fits'

def calculate_residuals(degree, n_frames=20):
    """Calculate residuals and save statistics for a given Legendre polynomial degree."""
    # Define paths for the current degree
    fit_cube_path = f'F:/legfit/{n_frames}_frames/fit_cube_leg_{n_frames}_frames{degree}deg_final_vst.fits'
    residuals_cube_path = f'F:/legfit/{n_frames}_frames/res_cube_leg_{n_frames}_frames{degree}deg_final_vst.fits'
    stat_table = f'F:/legfit/{n_frames}_frames/frame_statistics_leg_{n_frames}_frames{degree}deg.csv'

    # Load data
    y_cube_sliced = fits.getdata(y_cube_path)[1:n_frames]  # Load y_cube data
    y_cube = y_cube_sliced[:, 0, :, :]  # Take out the second dimension
    fit_cube = fits.getdata(fit_cube_path)  # Load fit cube data

    # Calculate residuals
    residuals_cube = y_cube - fit_cube
    fits.writeto(residuals_cube_path, residuals_cube, overwrite=True)  # Save residuals cube

    # Initialize lists to store statistics
    frame_num = []
    means = []
    rms_vals = []
    median_vals = []
    std_vals = []

    # Calculate statistics for each frame
    res = fits.getdata(residuals_cube_path)
    initial_frame_label = 1124972

    for i in range(res.shape[0]):
        data = res[i]
        means.append(np.mean(data))  # Calculate mean
        rms_vals.append(np.sqrt(np.mean(data ** 2)))  # Calculate RMS of residuals
        median_vals.append(np.median(data))  # Calculate median
        std_vals.append(np.std(data))  # Calculate std of residuals
        frame_num.append(initial_frame_label + i)  # Adjusted frame numbering

    # Save statistics to CSV
    table = pd.DataFrame({'Frame': frame_num, 'Mean': means, 'RMS': rms_vals, 'Median': median_vals, 'StdDev': std_vals})
    table.to_csv(stat_table, index=False)

    return residuals_cube

# Loop through degrees from 1 to 10
for degree in range(1, 11):
    calculate_residuals(degree)
