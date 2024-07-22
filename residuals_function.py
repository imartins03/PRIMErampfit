from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#definition of paths
super_bias_path = 'IRRC_calfiles\\super_biasC1.fits.ramp.20231012'
calFile = r'IRRC_calfiles\irrc_weights_C1.h5'
maskFile_path = r'IRRC_calfiles\C1_bad_ref_pix_mask.fits'
y_cube_path = r'D:\NLC\C1\y_cube_25.fits'
fit_cube_path = r'D:\NLC\C1\fit_cube.fits'
fit_coeff_path = r'D:\NLC\C1\fit_coeff.fits'
residuals_cube_path = r'D:\NLC\C1\residuals.fits'

def calculate_residuals():
    y_cube = fits.getdata(y_cube_path)  # Load y_cube data
    y_cube = y_cube[:, 0, :, :]  # Select the first element along the second axis
    fit_cube = fits.getdata(fit_cube_path)  # Load fit cube data
    residuals_cube = y_cube - fit_cube  # Calculate residuals
    fits.writeto(residuals_cube_path, residuals_cube, overwrite=True)  # Save residuals cube
    return residuals_cube

res = calculate_residuals()  # Calculate residuals

initial_frame_label = 1124972

#creating lists to append values for the table
frame_num = []
means = []
rms_vals = []
median_vals = []
std_vals = []
for i in range(res.shape[0]):
    data = res[i]
    means.append(np.mean(data))  # Calculate mean
    rms_vals.append(np.sqrt(np.mean(data ** 2)))  # Calculate RMS of residuals
    median_vals.append(np.median(data))  # Calculate median
    std_vals.append(np.std(data))  # Calculate std of residuals
    frame_num.append(initial_frame_label + i)  # Adjusted frame numbering

table = pd.DataFrame({'Mean': means, 'RMS': rms_vals, 'Median': median_vals, 'StdDev': std_vals})
table.to_csv(r'D:\NLC\C1\frame_statistics.csv', index=False)  # Save statistics to CSV

std = np.nanstd(res)
for i in range(res.shape[0]):
    residuals_frame = res[i]
    bins = np.arange(-2 * std, 2 * std, std / 20)  # Generate histogram bins
    hist = np.histogram(residuals_frame[np.isfinite(residuals_frame)], bins=bins)
    plt.bar(hist[1][:-1], hist[0], color='blue')  # Plot histogram
    plt.title(f'Histogram of Residuals for Frame {i + 1}')
    plt.xlabel('Residual Value')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.clf()
    plt.savefig(f'D:\\NLC\\C1\\hist_{i}.png')  # Save histogram plot
