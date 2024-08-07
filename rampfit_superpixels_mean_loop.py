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


def calculate_rms_square(data):
    return np.sqrt(np.mean(np.square(data))) ** 2


def generate_fit_cube(degrees, centers):
    y_cube = fits.getdata(y_cube_path)[1:n_frames]
    x, _, y, z = y_cube.shape
    y_cube = y_cube[:, 0, :, :]

    degree_list = []
    center_list = []
    slope_list = []

    time = np.arange(y_cube.shape[0], dtype=np.double)

    # Create a DataFrame to store average RMS for each superpixel, degree, and center
    avg_rms_data = []

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

            # Append results to lists
            degree_list.append(degree)
            center_list.append(np.flip(center))
            slope_list.append(slope)  # Store slope of polynomial fit

            # Create the center_str for file naming
            center_str = f"{center[1]}_{center[0]}"

            coeff_filename = os.path.join(data_directory,
                                          f"{fit_coeff_path}center_{center_str}_{degree}deg_{n_frames}frames_noframe1.fits")
            residuals_filename = os.path.join(data_directory,
                                              f"{residuals_base_path}center_{center_str}_{degree}deg_{n_frames}frames_noframe1.fits")

            fits.writeto(coeff_filename, coefficients, overwrite=True)
            residuals = y - fit
            fits.writeto(residuals_filename, residuals, overwrite=True)

            # Calculate RMS for each frame
            residuals_cube = fits.getdata(residuals_filename)
            rms_values = [np.sqrt(np.mean(np.square(frame))) for frame in residuals_cube]

            # Calculate squared RMS values and their average
            rms_square_values = [rms ** 2 for rms in rms_values]
            sum_squared_rms = np.sum(rms_square_values)
            rms_of_avg = np.sqrt(sum_squared_rms) / np.sqrt(len(rms_values) - 1)

            print(len(rms_values))

            avg_rms_data.append({
                'Degree': degree,
                'Center': np.flip(center),
                'RMS of the Average': rms_of_avg
            })

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

    # Create the main results DataFrame (with slope)
    results_df = pd.DataFrame({
        'Degree': degree_list,
        'Center': center_list,
        'Slope': slope_list
    })

    # Save the main results DataFrame to a CSV file
    csv_filename = os.path.join(data_directory, 'results_vs_degree_trial.csv')
    results_df.to_csv(csv_filename, index=True)
    print(f"Main results saved to {csv_filename}")

    # Create the average RMS DataFrame
    avg_rms_df = pd.DataFrame(avg_rms_data)
#f
    # Save the average RMS DataFrame to a CSV file
    avg_rms_csv_filename = os.path.join(data_directory, 'average_rms_vs_degree_and_center.csv')
    avg_rms_df.to_csv(avg_rms_csv_filename, index=True)
    print(f"Average RMS results saved to {avg_rms_csv_filename}")

    # Print both DataFrames
    print("Main Results DataFrame:")
    print(results_df)
    print("Average RMS DataFrame:")
    print(avg_rms_df)

    results_df.to_csv(csv_filename, index=True)
    print(f"Results saved to {csv_filename}")

    print(results_df)


# Parameters
centers = [(500, 2048), (2048, 2048), (3500, 2048)]  # Centers of the superpixels
degrees = np.linspace(1, 10, 10)
generate_fit_cube(degrees, centers)