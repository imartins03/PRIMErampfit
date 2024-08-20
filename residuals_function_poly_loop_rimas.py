import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.io import fits
import os


y_cube_path = r"G:\20240820.rimas.0010-0109.HK.mean.fits"

# fit_cube_path_template = r'D:\NLC\C1\dif_degrees_test\100_frames\fit_cube_poly_{degree}deg_{n_frames}frames_noframe1_weights.fits'
# fit_coeff_path_template = r'D:\NLC\C1\dif_degrees_test\100_frames\fit_coeff_poly_{degree}deg_{n_frames}frames_noframe1_weights.fits'
# residuals_cube_path_template = r'F:\leftover_C1_dif_degrees_test_rampfit\239_frames\residuals_poly_{degree}deg_239frames_noframe1.fits'

#239 frame paths
fit_cube_path_template = r"G:\fit_cube_poly_{degree}deg_{n_frames}frames.fits"
fit_coeff_path_template = r"G:\fit_coeff_poly_{degree}deg_{n_frames}frames.fits"
residuals_cube_path_template = r"G:\residuals_poly_{degree}deg_{n_frames}frames.fits"

frame_statistics_table_template = r"G:\frame_statistics_poly_{degree}deg_{n_frames}frames.csv"
slope_statistics_table_path = r"G:\slope_statistics.csv"
rms_statistics_table_path = r"G:\rms_statistics.csv"

n_frames = 22

def calculate_residuals(degree):
    y_cube = fits.getdata(y_cube_path)

    fit_cube_path = fit_cube_path_template.format(degree=degree,n_frames=n_frames)
    residuals_cube_path = residuals_cube_path_template.format(degree=degree,n_frames=n_frames)

    fit_cube = fits.getdata(fit_cube_path)
    fit_coeff = fits.getdata(fit_coeff_path_template.format(degree=degree,n_frames=n_frames))

    residuals_cube = y_cube - fit_cube  # Calculate residuals
    fits.writeto(residuals_cube_path, residuals_cube, overwrite=True)  # Save residuals cube
    return residuals_cube, fit_coeff, fit_cube

def compute_statistics(residuals_cube, fit_coeff, initial_frame_label):
    means = []
    rms_vals = []
    median_vals = []
    std_vals = []
    slopes = []
    frame_num = []

    slope_vals = fit_coeff[-2]  # Extract the slope layer

    for i in range(residuals_cube.shape[0]):
        data = residuals_cube[i]
        means.append(np.mean(data))
        rms_vals.append(np.sqrt(np.mean(data ** 2)))
        median_vals.append(np.median(data))
        std_vals.append(np.std(data))

        slopes.append(np.mean(slope_vals))  # Store mean slope value for the current frame
        frame_num.append(initial_frame_label + i)

    # Save frame statistics
    frame_stats_df = pd.DataFrame({
        'FrameNumber': frame_num,
        'Mean': means,
        'RMS': rms_vals,
        'Median': median_vals,
        'StdDev': std_vals
    })
    frame_stats_df.to_csv(frame_statistics_table_template.format(degree=degree,n_frames=n_frames, index=False))

    # Aggregate statistics
    rms_of_avg = np.sqrt(np.sum(np.array(rms_vals) ** 2))/np.sqrt(len(rms_vals)-1)
    avg_slope = np.mean(slopes)
    slope_sem = np.std(slope_vals) / np.sqrt(4088 * 4088)  # Updated SEM calculation

    return rms_of_avg, avg_slope, slope_sem

def save_statistics(degree, rms_of_avg, avg_slope, slope_sem):
    # Save RMS statistics
    rms_stats_df = pd.DataFrame({
        'DegreeOfFit': [degree],
        'RMSofAverage': [rms_of_avg]
    })
    rms_stats_df.to_csv(rms_statistics_table_path, mode='a', header=not os.path.exists(rms_statistics_table_path), index=False)
#
    # Save slope statistics
    slope_stats_df = pd.DataFrame({
        'DegreeOfFit': [degree],
        'SlopeOfFit': [avg_slope],
        'SlopeSEM': [slope_sem]
    })
    slope_stats_df.to_csv(slope_statistics_table_path, mode='a', header=not os.path.exists(slope_statistics_table_path), index=False)

def plot_statistics():
    df_rms = pd.read_csv(rms_statistics_table_path)
    df_slope = pd.read_csv(slope_statistics_table_path)

    plt.figure()
    plt.plot(df_rms['DegreeOfFit'], df_rms['RMSofAverage'], marker='o', linestyle='-')
    plt.xlabel('Degree of Fit')
    plt.ylabel('Average RMS')
    plt.title('Average RMS vs Degree of Fit')
    plt.grid(True)
    plt.savefig(r"F:\leftover_C1_dif_degrees_test_rampfit\239_frames\weighted\average_rms_vs_degree.png")
    plt.close()

    plt.figure()
    plt.plot(df_slope['DegreeOfFit'], df_slope['SlopeOfFit'], marker='o', linestyle='-')
    plt.xlabel('Degree of Fit')
    plt.ylabel('Average Slope')
    plt.title('Average Slope vs Degree of Fit')
    plt.grid(True)

    plt.savefig(r"F:\leftover_C1_dif_degrees_test_rampfit\239_frames\weighted\average_slope_vs_degree.png")
    plt.close()

    plt.figure()
    plt.plot(df_slope['DegreeOfFit'], df_slope['SlopeSEM'], marker='o', linestyle='-')
    plt.xlabel('Degree of Fit')
    plt.ylabel('Slope SEM')
    plt.title('Slope SEM vs Degree of Fit')
    plt.grid(True)

    plt.savefig(r"F:\leftover_C1_dif_degrees_test_rampfit\239_frames\weighted\slope_sem_vs_degree.png")
    plt.close()

initial_frame_label = 1124973

total_degrees = 10
for degree in range(1, total_degrees + 1):
    print(f"Processing degree {degree}/{total_degrees} ({(degree / total_degrees) * 100:.2f}%)")
    residuals_cube, fit_coeff, fit_cube = calculate_residuals(degree)
    rms_of_avg, avg_slope, slope_sem = compute_statistics(residuals_cube, fit_coeff, initial_frame_label)
    save_statistics(degree, rms_of_avg, avg_slope, slope_sem)

plot_statistics()
##