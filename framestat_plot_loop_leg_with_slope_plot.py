import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
import os

# Paths
y_cube_path = r'D:\NLC\C1\y_cube_500.fits'
fit_cube_path_template = r"F:\legfit\239_frames_unweighted\final_coefficient_and_fit_images\fit_cube_leg_239frames_{degree}deg_final_vst_unweighted.fits"
residuals_cube_path_template = r"F:\legfit\239_frames_unweighted\res_cube_leg_239frames_{degree}deg_final_vst_unweighted.fits"
fit_coeff_path_template = r"F:\legfit\239_frames_unweighted\final_coefficient_and_fit_images\coefficients_leg_239frames_{degree}deg_final_vst_unweighted.fits"
stat_table_template = r'F:\legfit\239_frames_unweighted\res_stats\frame_statistics_leg_{degree}deg_239frames_noframe1.csv'
table_path = r'F:\legfit\239_frames_unweighted\res_stats\legfit_res_rms_avg_unweighted.csv'

degree_range = np.arange(1, 11)
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
    frame_num = []
    slopes = []

    slope_vals = fit_coeff[-2]  # Assuming this is the layer with slope data

    for i in range(residuals_cube.shape[0]):
        data = residuals_cube[i]
        means.append(np.mean(data))
        rms_vals.append(np.sqrt(np.mean(data ** 2)))
        median_vals.append(np.median(data))
        std_vals.append(np.std(data))
        slopes.append(np.mean(slope_vals[i]))  # Mean slope for this frame
        frame_num.append(initial_frame_label + i)

    frame_stats_df = pd.DataFrame({
        'FrameNumber': frame_num,
        'Mean': means,
        'RMS': rms_vals,
        'Median': median_vals,
        'StdDev': std_vals,
    })

    frame_stats_df.to_csv(stat_table_template.format(degree=degree), index=False)

    total_rms_square_sum = np.sum(np.array(rms_vals) ** 2)
    length_of_data = len(rms_vals)
    divisor = np.sqrt(length_of_data - 1)
    rms_of_avg = np.sqrt(total_rms_square_sum) / divisor

    mean_slope = np.mean(slopes)  # Mean slope across all frames

    return rms_of_avg, mean_slope

def save_statistics(degree, rms_of_avg, mean_slope):
    rms_slope_df = pd.DataFrame({
        'DegreeOfFit': [degree],
        'RMSofAverage': [rms_of_avg],
        'SlopeofFit': [mean_slope]
    })

    rms_slope_df.to_csv(table_path, mode='a', header=not os.path.exists(table_path), index=False)

def plot_statistics():
    df_rms_slope = pd.read_csv(table_path)

    plt.figure(1)
    plt.plot(df_rms_slope['DegreeOfFit'], df_rms_slope['RMSofAverage'], marker='o', linestyle='-')
    plt.xlabel('Degree of Fit')
    plt.ylabel('Average RMS')
    plt.title('Average RMS vs Degree of Fit')
    plt.grid(True)
    plt.savefig(r'F:\legfit\239_frames_unweighted\res_stats\average_rms_vs_degree_legfit_239frames_noframe1_unweighted.png')

    plt.figure(2)
    plt.plot(df_rms_slope['DegreeOfFit'], df_rms_slope['SlopeofFit'], marker='o', linestyle='-')
    plt.xlabel('Degree of Fit')
    plt.ylabel('Slope of Fit')
    plt.title('Slope of Fit vs Degree of Fit')
    plt.grid(True)
    plt.savefig(r'F:\legfit\239_frames_unweighted\res_stats\slope_of_the_fit_vs_degree_legfit_239_frames_noframe1_unweighted.png')

initial_frame_label = 1124973  # Start frame number

for degree in degree_range:
    print(f"Processing degree {degree}...")
    residuals_cube, fit_cube, fit_coeff = load_data(degree)
    rms_of_avg, mean_slope = compute_statistics(residuals_cube, fit_coeff, initial_frame_label)
    save_statistics(degree, rms_of_avg, mean_slope)

plot_statistics()  # Plot and save the statistics plots
