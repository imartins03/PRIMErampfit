from astropy.io import fits
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

n_frames = 239

# Definition of paths
super_bias_path = 'IRRC_calfiles\\super_biasC1.fits.ramp.20231012'
calFile = r'IRRC_calfiles\irrc_weights_C1.h5'
maskFile_path = r'IRRC_calfiles\C1_bad_ref_pix_mask.fits'
y_cube_path = r'D:\NLC\C1\y_cube_500.fits'
fit_cube_path_template = r'F:\leftover_C1_dif_degrees_test_rampfit\239_frames\fit_cube_poly_{degree}deg_239frames_noframe1.fits'
fit_coeff_path_template = r'F:\leftover_C1_dif_degrees_test_rampfit\239_frames\fit_coeff_poly_{degree}deg_239frames_noframe1.fits'
residuals_cube_path_template = r'F:\leftover_C1_dif_degrees_test_rampfit\239_frames\residuals_poly_{degree}deg_239frames_noframe1.fits'
stat_table_template = r'F:\leftover_C1_dif_degrees_test_rampfit\239_frames\frame_statistics_poly_{degree}deg_239frames_noframe1.csv'

def calculate_residuals(degree):
    y_cube = fits.getdata(y_cube_path)[1:n_frames]  # Load y_cube data
    y_cube = y_cube[:,0, :, :]  # Slice out frames 2 to 100
    fit_cube_path = fit_cube_path_template.format(degree=degree)
    residuals_cube_path = residuals_cube_path_template.format(degree=degree)

    fit_cube = fits.getdata(fit_cube_path)[:n_frames]  # Load fit cube data
    fit_coeff = fits.getdata(fit_coeff_path_template.format(degree=degree))

    residuals_cube = y_cube - fit_cube  # Calculate residuals
    fits.writeto(residuals_cube_path, residuals_cube, overwrite=True)  # Save residuals cube
    return residuals_cube, fit_coeff

def compute_statistics(residuals_cube, fit_coeff, initial_frame_label):
    means = []
    rms_vals = []
    median_vals = []
    std_vals = []
    slopes = []  # List to store slope values
    frame_num = []

    # Assuming fit_coeff shape is (2, 4088, 4088) for a first-degree polynomial fit
    slope_vals = fit_coeff[1]  # Extract the slope layer

    for i in range(residuals_cube.shape[0]):
        data = residuals_cube[i]
        means.append(np.mean(data))  # Calculate mean
        rms_vals.append(np.sqrt(np.mean(data ** 2)))  # Calculate RMS of residuals
        median_vals.append(np.median(data))  # Calculate median
        std_vals.append(np.std(data))  # Calculate std of residuals

        # For each frame, the slope can be extracted from the slope_coefficients
        # Assuming you need a single slope value per frame
        # Here, we are assuming the slope is consistent across the frame.
        slopes.append((slope_vals))  # Example: average slope across the entire frame

        frame_num.append(initial_frame_label + i)  # Adjusted frame numbering

    table = pd.DataFrame({
        'FrameNumber': frame_num,
        'Mean': means,
        'RMS': rms_vals,
        'Median': median_vals,
        'StdDev': std_vals,
        'Slope': slopes,  # Add Slope column
    })

    return table

def save_statistics(df, degree):
    total_rms_square_sum = np.sum((df['RMS'])**2)
    length_of_data = len(df)
    divisor = np.sqrt(length_of_data - 1)
    df['TotalRMS Square Sum'] = total_rms_square_sum
    df['RMS of the Average'] = df['TotalRMS'] / divisor
    stat_table = stat_table_template.format(degree=degree)
    df.to_csv(stat_table, index=False)

# Process for degrees from 1 to 10
initial_frame_label = 1124973 #start one later since first frame was cut out

for degree in range(1, 11):
    print(f"Processing degree {degree}")

    residuals_cube, fit_coeff = calculate_residuals(degree)  # Calculate residuals
    statistics_df = compute_statistics(residuals_cube, fit_coeff, initial_frame_label)  # Compute statistics
    save_statistics(statistics_df, degree)  # Save statistics to CSV

    # Plot histograms for each frame
    # for i in range(residuals_cube.shape[0]):
    #     residuals_frame = residuals_cube[i]
    #     mean = statistics_df.loc[i, 'Mean']
    #     std = statistics_df.loc[i, 'StdDev']
    #     bins = np.arange(-3 * std + mean, mean + 3 * std, std / 20)  # Generate histogram bins
    #     hist = np.histogram(residuals_frame[np.isfinite(residuals_frame)], bins=bins)
    #     plt.bar(hist[1][:-1], hist[0], color='blue', width=std / 20)  # Plot histogram
    #     plt.title(f'Histogram of Residuals for Frame {i + 1}')
    #     plt.xlabel('Residual Value')
    #     plt.ylabel('Frequency')
    #     plt.xlim((bins.min(), bins.max()))
    #     plt.grid(True)
    #     plt.savefig(f'D:\\NLC\\C1\\hist_degree_{degree}_frame_{i}_noframe1.png')  # Save histogram plot
    #     if i % 10 == 0:
    #         plt.show()
    #     plt.clf()
