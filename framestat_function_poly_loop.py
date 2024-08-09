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
rms_slope_statistics_table_path = r'F:\leftover_C1_dif_degrees_test_rampfit\rms_slope_statistics_attempt.csv'

n_frames = 239

def load_data(degree):
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

    slope_vals = fit_coeff[-2]  # Extract the slope layer

    # for i in range(residuals_cube.shape[0]):
    #     data = residuals_cube[i]
    #     means.append(np.mean(data))
    #     rms_vals.append(np.sqrt(np.mean(data ** 2)))
    #     median_vals.append(np.median(data))
    #     std_vals.append(np.std(data))
    #
    #     slopes.append(np.mean(slope_vals))  # Store mean slope value for the current frame
    #     frame_num.append(initial_frame_label + i)
    #
    # # Save frame statistics
    # frame_stats_df = pd.DataFrame({
    #     'FrameNumber': frame_num,
    #     'Mean': means,
    #     'RMS': rms_vals,
    #     'Median': median_vals,
    #     'StdDev': std_vals
    # })
    # frame_stats_df.to_csv(stat_table_template.format(degree=degree), index=False)

    # Aggregate statistics
    rms_of_avg = np.sqrt(np.mean(np.array(rms_vals) ** 2))
    avg_slope = np.mean(np.mean(slope_vals))
    slope_sem = sem(slopes)  # Calculate SEM for the slopes

    return rms_of_avg, avg_slope, slope_sem

def save_rms_of_average_and_slope_statistics(degree, rms_of_avg, avg_slope, slope_sem):
    file_exists = os.path.exists(rms_slope_statistics_table_path)
    rms_slope_df = pd.DataFrame({
        'DegreeOfFit': [degree],
        'RMSofAverage': [rms_of_avg],
        'SlopeOfFit': [avg_slope],
        'SlopeSEM': [slope_sem]
    })
    rms_slope_df.to_csv(rms_slope_statistics_table_path, mode='a', header=not file_exists, index=False)

def plot_statistics():
    df_rms_slope = pd.read_csv(rms_slope_statistics_table_path)

    plt.figure()
    plt.plot(df_rms_slope['DegreeOfFit'], df_rms_slope['RMSofAverage'], marker='o', linestyle='-')
    plt.xlabel('Degree of Fit')
    plt.ylabel('Average RMS')
    plt.title('Average RMS vs Degree of Fit')
    plt.grid(True)
    plt.savefig('average_rms_vs_degree.png')
    plt.close()

    plt.figure()
    plt.plot(df_rms_slope['DegreeOfFit'], df_rms_slope['SlopeOfFit'], marker='o', linestyle='-')
    plt.xlabel('Degree of Fit')
    plt.ylabel('Average Slope')
    plt.title('Average Slope vs Degree of Fit')
    plt.grid(True)
    plt.savefig('average_slope_vs_degree.png')
    plt.close()

    plt.figure()
    plt.plot(df_rms_slope['DegreeOfFit'], df_rms_slope['SlopeSEM'], marker='o', linestyle='-')
    plt.xlabel('Degree of Fit')
    plt.ylabel('Slope SEM')
    plt.title('Slope SEM vs Degree of Fit')
    plt.grid(True)
    plt.savefig('slope_sem_vs_degree.png')
    plt.close()

initial_frame_label = 1124973

total_degrees = 10
for degree in range(1, total_degrees + 1):
    print(f"Processing degree {degree}/{total_degrees} ({(degree / total_degrees) * 100:.2f}%)")
    residuals_cube, fit_cube, fit_coeff = load_data(degree)
    rms_of_avg, avg_slope, slope_sem = compute_statistics(residuals_cube, fit_coeff, initial_frame_label)
    save_rms_of_average_and_slope_statistics(degree, rms_of_avg, avg_slope, slope_sem)

plot_statistics()
