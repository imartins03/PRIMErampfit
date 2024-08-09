import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
import os
from scipy.stats import sem

# Paths
y_cube_path = r'D:\NLC\C1\y_cube_500.fits'
fit_cube_path_template = r'F:\leftover_C1_dif_degrees_test_rampfit\239_frames\fit_cube_poly_{degree}deg_239frames_noframe1.fits'
residuals_cube_path_template = r'F:\leftover_C1_dif_degrees_test_rampfit\239_frames\residuals_poly_{degree}deg_239frames_noframe1.fits'
fit_coeff_path_template = r'F:\leftover_C1_dif_degrees_test_rampfit\239_frames\fit_coeff_poly_{degree}deg_239frames_noframe1.fits'
stat_table_template = r'F:\leftover_C1_dif_degrees_test_rampfit\239_frames\frame_statistics_poly_{degree}deg_239frames_noframe1_attempt.csv'
rms_slope_statistics_table_path = r'F:\leftover_C1_dif_degrees_test_rampfit\rms_slope_statistics_attempt_new2.csv'

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
    slopes = []
    frame_num = []

    # fit_coeff shape is (2, 4088, 4088) for a first-degree polynomial fit
    slope_vals = fit_coeff[-2]  # Extract the slope layer

    # for i in range(residuals_cube.shape[0]):
    #     data = residuals_cube[i]
    #     means.append(np.mean(data))  # Calculate mean
    #     rms_vals.append(np.sqrt(np.mean(data ** 2)))  # Calculate RMS of residuals
    #     median_vals.append(np.median(data))  # Calculate median
    #     std_vals.append(np.std(data))  # Calculate std of residuals
    #
    #     slopes.append(np.mean(slope_vals))
    #     frame_num.append(initial_frame_label + i)
    #
    # # Save statistics for frames
    # frame_stats_df = pd.DataFrame({
    #     'FrameNumber': frame_num,
    #     'Mean': means,
    #     'RMS': rms_vals,
    #     'Median': median_vals,
    #     'StdDev': std_vals
    # })
    #
    # # Save statistics to CSV, overwriting the file each time
    # frame_stats_df.to_csv(stat_table_template.format(degree=degree), index=False)

    # Compute RMS of average
    total_rms_square_sum = np.sum(np.array(rms_vals) ** 2)
    length_of_data = n_frames
    divisor = np.sqrt(length_of_data - 1)
    rms_of_avg = np.sqrt(total_rms_square_sum) / divisor
    avg_slope = np.mean(slope_vals)  # Compute the mean slope
    stddev = np.std(slope_vals)

    # Calculate the standard error of the mean (SEM) for slopes
    # slope_sem = sem(slope_vals)
    slope_sem = np.std(slope_vals)/np.sqrt(4088*4088)

    # return rms_of_avg, avg_slope, slope_sem,std
    return avg_slope, slope_sem, stddev
#
def save_rms_of_average_and_slope_statistics(degree, avg_slope, slope_sem, stddev):
    rms_slope_df = pd.DataFrame({
        'DegreeOfFit': [degree],
        # 'RMSofAverage': [rms_of_avg],
        'SlopeOfFit': [avg_slope],
        'SlopeSEM': [slope_sem],
        'stdev': [stddev]
    })

    # Save RMS and slope to CSV, appending to the file each time
    rms_slope_df.to_csv(rms_slope_statistics_table_path, mode='a', header=not os.path.exists(rms_slope_statistics_table_path), index=False)

def plot_statistics():
    # Load the statistics data for plotting
    df_rms_slope = pd.read_csv(rms_slope_statistics_table_path)

    # # Plot average RMS vs degree of fit
    # plt.figure()
    # plt.plot(df_rms_slope['DegreeOfFit'], df_rms_slope['RMSofAverage'], marker='o', linestyle='-')
    # plt.xlabel('Degree of Fit')
    # plt.ylabel('Average RMS')
    # plt.title('Average RMS vs Degree of Fit')
    # plt.grid(True)
    # plt.savefig('average_rms_vs_degree.png')
    # plt.close()

    # Plot slope of the fit vs degree of fit
    plt.figure()
    plt.plot(df_rms_slope['DegreeOfFit'], df_rms_slope['SlopeOfFit'], marker='o', linestyle='-')
    plt.xlabel('Degree of Fit')
    plt.ylabel('Average Slope')
    plt.title('Average Slope vs Degree of Fit')
    plt.grid(True)
    plt.savefig('average_slope_vs_degree.png')
    plt.close()

    # Plot SEM of slope vs degree of fit
    plt.figure()
    plt.plot(df_rms_slope['DegreeOfFit'], df_rms_slope['SlopeSEM'], marker='o', linestyle='-')
    plt.xlabel('Degree of Fit')
    plt.ylabel('Slope SEM')
    plt.title('Slope SEM vs Degree of Fit')
    plt.grid(True)
    plt.savefig('slope_sem_vs_degree.png')
    plt.close()

initial_frame_label = 1124973  # start one later since the first frame was cut out

total_degrees = 10
for degree in range(1, total_degrees + 1):
    print(f"Processing degree {degree}/{total_degrees} ({(degree / total_degrees) * 100:.2f}%)")

    residuals_cube, fit_cube, fit_coeff = load_data(degree)  # Load data
    avg_slope, slope_sem, stddev = compute_statistics(residuals_cube, fit_coeff, initial_frame_label)  # Compute statistics

    save_rms_of_average_and_slope_statistics(degree, avg_slope, slope_sem, stddev)

plot_statistics()  # Plot and save the statistics plots