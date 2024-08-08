from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import matplotlib.patches as patches

# Paths
super_bias_path = 'IRRC_calfiles\\super_biasC1.fits.ramp.20231012'
maskFile_path = r'IRRC_calfiles\C1_bad_ref_pix_mask.fits'
y_cube_path = r'D:\NLC\C1\y_cube_500.fits'
data_directory = r'D:\NLC\C1\superpix\mean_superpixel'
fit_cube_base_path = r'fit_cube_avg_'
fit_coeff_path = r'fit_coeff_avg_'
residuals_base_path = r'residuals_avg_'

# Load data
super_bias = fits.getdata(super_bias_path)  # Load super bias data
maskFile = fits.getdata(maskFile_path)  # Load mask file data
mask = maskFile > 0  # Create mask for bad pixels
supercpy = super_bias.copy()  # Create a copy of the super bias
supercpy[mask[:, :4096]] = 0  # Apply mask to the super bias copy
n_frames = 239


def generate_fit_cube(degrees, centers):
    # Define the base paths for saving files
    fit_coeff_path = r'fit_coeff_avg_'
    residuals_base_path = r'residuals_avg_'
    data_directory = r'D:\NLC\C1\superpix\mean_superpixel'

    # Load data
    y_cube = fits.getdata(y_cube_path)[1:n_frames]
    x, _, y, z = y_cube.shape
    y_cube = y_cube[:, 0, :, :]

    # Initialize lists and dictionaries for tracking data
    combined_results_data = []
    slopes_by_center = {center: [] for center in centers}
    time = np.arange(y_cube.shape[0], dtype=np.double)

    for degree in degrees:
        for center in centers:
            half_size = 128
            y_start = center[0] - half_size
            y_end = center[0] + half_size
            x_start = center[1] - half_size
            x_end = center[1] + half_size
            region = (y_start, y_end, x_start, x_end)

            y = y_cube[:, region[0]:region[1], region[2]:region[3]]
            y = y.reshape(x, -1)
            y = np.mean(y, axis=1)

            # Polynomial fitting
            coefficients = np.polyfit(time, y, degree)
            fit = np.polyval(coefficients, time)

            slope = coefficients[-2]

            # Calculate RMS for each frame
            residuals = y - fit
            rms_values = [np.sqrt(np.mean(np.square(frame))) for frame in residuals]

            # Calculate squared RMS values and their average
            rms_square_values = [rms ** 2 for rms in rms_values]
            sum_squared_rms = np.sum(rms_square_values)
            rms_of_avg = np.sqrt(sum_squared_rms) / np.sqrt(len(rms_values) - 1)

            # Append results to the combined results list
            combined_results_data.append({
                'Degree': degree,
                'Center': np.flip(center),
                'RMS of the Average': rms_of_avg,
                'Slope': slope
            })

            # Track slopes for standard deviation calculation
            slopes_by_center[tuple(center)].append(slope)

            # Create the center_str for file naming
            center_str = f"{center[1]}_{center[0]}"

            # Define file paths for coefficients and residuals
            coeff_filename = os.path.join(data_directory,
                                          f"{fit_coeff_path}center_{center_str}_{degree}deg_{n_frames}frames_noframe1.fits")
            residuals_filename = os.path.join(data_directory,
                                              f"{residuals_base_path}center_{center_str}_{degree}deg_{n_frames}frames_noframe1.fits")

            # Save coefficients and residuals as FITS files
            fits.writeto(coeff_filename, coefficients, overwrite=True)
            fits.writeto(residuals_filename, residuals, overwrite=True)

        # Plot the residuals
        plt.figure()
        for center in centers:
            center_str = f"{center[1]}_{center[0]}"
            residuals_filename = os.path.join(data_directory,
                                              f"{residuals_base_path}center_{center_str}_{degree}deg_{n_frames}frames_noframe1.fits")
            residuals = fits.getdata(residuals_filename)
            plt.plot(time, residuals, label=f'Center {np.flip(center)}')

        plt.title(f'Residuals vs. Frame Number for Degree {degree}, Avg Superpixel, no frame 1')
        plt.xlabel('Frame Number')
        plt.ylabel('Residuals')
        plt.legend()
        plt.grid(True)

        plot_filename = os.path.join(data_directory,
                                     f'median_avg_superpix_plot_poly_{degree}deg_{n_frames}frames_noframe1.png')
        plt.savefig(plot_filename)
        plt.show()
        plt.close()

    # Create the combined results DataFrame
    combined_results_df = pd.DataFrame(combined_results_data)

    # Save the combined results DataFrame to a CSV file
    combined_results_csv_filename = os.path.join(data_directory, 'combined_results.csv')
    combined_results_df.to_csv(combined_results_csv_filename, index=False)
    print(f"Combined results saved to {combined_results_csv_filename}")

    # Calculate the standard deviation of slopes for each center
    stddev_slope_data = []
    for center, slopes in slopes_by_center.items():
        stddev_slope = np.std(slopes)

        stddev_slope_data.append({
            'Center': np.flip(center),
            'Standard Deviation of Slope': stddev_slope
        })

    # Create the standard deviation DataFrame
    stddev_slope_df = pd.DataFrame(stddev_slope_data)

    # Save the standard deviation DataFrame to a CSV file
    stddev_slope_csv_filename = os.path.join(data_directory, 'stddev_slope_vs_center.csv')
    stddev_slope_df.to_csv(stddev_slope_csv_filename, index=False)
    print(f"Standard deviation of slopes results saved to {stddev_slope_csv_filename}")

    # Print both DataFrames
    print("Combined Results DataFrame:")
    print(combined_results_df)
    print("Standard Deviation of Slopes DataFrame:")
    print(stddev_slope_df)


# Parameters
centers = [(500, 2048), (2048, 2048), (3500, 2048)]  # Centers of the superpixels
degrees = np.linspace(1, 10, 10)  # Degrees of polynomial fit
generate_fit_cube(degrees, centers)
