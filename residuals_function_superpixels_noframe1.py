from astropy.io import fits
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Paths
y_cube_path = r'D:\NLC\C1\y_cube_4.fits'
fit_cube_base_path = r'D:\NLC\C1\fit_cube_'
residuals_cube_base_path = r'D:\NLC\C1\residuals_'
stat_table_base_path = r'D:\NLC\C1\frame_statistics_'

# Superpixel centers and parameters
centers = [(2048, 2048), (3072, 2048), (1024, 2048)]  # Centers of the superpixels
degrees = 1
n_frames = 100
initial_frame_label = 1124973

def calculate_residuals_for_superpixel(center, degrees):
    center_str = f"{center[1]}_{center[0]}"
    fit_cube_path = fit_cube_base_path + f'center_{center_str}_{degrees}deg_noframe1.fits'
    residuals_cube_path = residuals_cube_base_path + f'center_{center_str}_{degrees}deg_noframe1.fits'
    stat_table_path = stat_table_base_path + f'center_{center_str}_{degrees}deg_noframe1.csv'

    # "D:\NLC\C1\fit_cube_center_2048_2048_1deg_noframe1.fits"

    # Load y_cube and fit_cube
    y_cube = fits.getdata(y_cube_path)[1:n_frames]  # Load y_cube data
    y_cube = y_cube[:, 0, :, :]  # Remove the second dimension


    # Define region
    half_size = 128  # Size for superpixel
    y_start = center[0] - half_size
    y_end = center[0] + half_size
    x_start = center[1] - half_size
    x_end = center[1] + half_size

    # Slice y_cube and fit_cube to superpixel region
    y_cube = y_cube[:, y_start:y_end, x_start:x_end]
    fit_cube = fits.getdata(fit_cube_path)[:n_frames]
    print(np.shape(y_cube), np.shape(fit_cube))

    # Calculate residuals
    residuals_cube = y_cube - fit_cube
    fits.writeto(residuals_cube_path, residuals_cube, overwrite=True)  # Save residuals cube

    # Create statistics table
    frame_num = []
    means = []
    rms_vals = []
    median_vals = []
    std_vals = []
    error = []  # Mean(rms)/sqrt(n_frames)

    for i in range(residuals_cube.shape[0]):
        data = residuals_cube[i]
        means.append(np.mean(data))  # Calculate mean
        rms_vals.append(np.sqrt(np.mean(data ** 2)))  # Calculate RMS of residuals
        median_vals.append(np.median(data))  # Calculate median
        std_vals.append(np.std(data))  # Calculate std of residuals
        error.append(np.mean(rms_vals) / np.sqrt(len(data)))
        frame_num.append(initial_frame_label + i)  # Adjusted frame numbering

    # Save statistics to CSV
    table = pd.DataFrame({
        'Frame': frame_num,
        'Mean': means,
        'RMS': rms_vals,
        'Median': median_vals,
        'StdDev': std_vals,
        'Error': error
    })
    table.to_csv(stat_table_path, index=False)

    return residuals_cube

# Process each superpixel and generate residuals and statistics
for center in centers:
    residuals_cube = calculate_residuals_for_superpixel(center, degrees)
    print(f'Processed superpixel at center {center}')
