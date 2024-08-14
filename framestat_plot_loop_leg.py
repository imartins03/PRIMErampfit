import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
import os

# Paths
y_cube_path = r'D:\NLC\C1\y_cube_500.fits'
fit_cube_path_template = r"F:\legfit\239_frames_unweighted\fit_cube_leg_239frames_{degree}deg_final_vst_unweighted.fits"
residuals_cube_path_template = r"F:\legfit\239_frames_unweighted\res_cube_leg_239frames_{degree}deg_final_vst_unweighted.fits"
fit_coeff_path_template = r"F:\legfit\239_frames_unweighted\coefficients_leg_239frames_{degree}deg_final_vst_unweighted.fits"
stat_table_template = r'"F:\legfit\239_frames_unweighted\res_stats\frame_statistics_leg_{degree}deg_239frames_noframe1.csv'
# Updated path for the second table
table_path = r'F:\leftover_C1_dif_degrees_test_rampfit\239_frames\frame_statistics.csv'

degree = np.arange(1,11)
n_frames = 239

def load_data(degree):
    # Load data from files
    residuals_cube = fits.getdata(residuals_cube_path_template.format(degree=degree))
    fit_cube = fits.getdata(fit_cube_path_template.format(degree=degree))[:n_frames]
    fit_coeff = fits.getdata(fit_coeff_path_template.format(degree=degree))

    return residuals_cube, fit_cube, fit_coeff

def compute_statistics(residuals_cube, fit_coeff, initial_frame_label):
    means = []
    rms_vals = []
    median_vals = []
    std_vals = []
    variances = []
    slopes = []
    frame_num = []

    # fit_coeff shape is (2, 4088, 4088) for a first-degree polynomial fit


    for i in range(residuals_cube.shape[0]):
        data = fits.getdata(residuals_cube[i])
        means.append(np.mean(data))  # Calculate mean
        rms_vals.append(np.sqrt(np.mean(data ** 2)))  # Calculate RMS of residuals
        median_vals.append(np.median(data))  # Calculate median
        std_vals.append(np.std(data))  # Calculate std of residuals
        frame_num.append(initial_frame_label + i)  # Adjusted frame numbering

    # Save statistics for frames
    frame_stats_df = pd.DataFrame({
        'FrameNumber': frame_num,
        'Mean': means,
        'RMS': rms_vals,
        'Median': median_vals,
        'StdDev': std_vals,
    })

    # Save statistics to CSV, overwriting the file each time
    frame_stats_df.to_csv(stat_table_template.format(degree=degree), index=False)

    # Compute RMS of average and mean slope
    total_rms_square_sum = np.sum(np.array(rms_vals) ** 2)
    length_of_data = len(rms_vals)
    divisor = np.sqrt(length_of_data - 1)
    rms_of_avg = np.sqrt(total_rms_square_sum) / divisor

    return rms_of_avg

def save_rms_of_average_and_slope_statistics(degree, rms_of_avg):
    rms_slope_df = pd.DataFrame({
        'DegreeOfFit': [degree],
        'RMSofAverage': [rms_of_avg],
    })

    # Save RMS and slope to CSV, appending to the file each time
    rms_slope_df.to_csv(table_path, mode='a', header=not os.path.exists(table_path), index=False)

def plot_statistics():
    # Load the statistics data for plotting
    df_rms_slope = pd.read_csv(table_path)

    # Plot average RMS vs degree of fit
    plt.figure()
    plt.plot(df_rms_slope['DegreeOfFit'], df_rms_slope['RMSofAverage'], marker='o', linestyle='-')
    plt.xlabel('Degree of Fit')
    plt.ylabel('Average RMS')
    plt.title('Average RMS vs Degree of Fit')
    plt.grid(True)

    plt.savefig(r'F:\leftover_C1_dif_degrees_test_rampfit\average_rms_vs_degree.png')
    plt.show()


    # Plot slope of the fit vs degree of fit
    # plt.figure()
    # plt.plot(df_rms_slope['DegreeOfFit'], df_rms_slope['SlopeOfFit'], marker='o', linestyle='-')
    # plt.xlabel('Degree of Fit')
    # plt.ylabel('Average Slope')
    # plt.title('Average Slope vs Degree of Fit')
    # plt.grid(True)
    #
    # plt.savefig(r'F:\leftover_C1_dif_degrees_test_rampfit\average_slope_vs_degree.png')
    # plt.show()

initial_frame_label = 1124973  # start one later since the first frame was cut out

total_degrees = 10

load_data(degree=degree)
compute_statistics(residuals_cube, fit_coeff, initial_frame_label)
plot_statistics()  # Plot and save the statistics plots
