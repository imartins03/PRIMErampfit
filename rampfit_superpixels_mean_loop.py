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


def calculate_rms(data):
    return np.sqrt(np.mean(np.square(data)))


def generate_fit_cube(degrees, centers):
    y_cube = fits.getdata(y_cube_path)[1:n_frames]
    x, _, y, z = y_cube.shape
    y_cube = y_cube[:, 0, :, :]

    degree_list = []
    center_list = []
    slope_list = []
    rms_list = []

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

            # Polynomial fitting and residuals
            coefficients = np.polyfit(time, y, degree)
            fit = np.polyval(coefficients, time)
            residuals = y - fit

            slope = (residuals[-1] - residuals[0]) / (time[-1] - time[0])

            # RMS of the residuals
            rms = calculate_rms(residuals)

            # Append results to lists
            degree_list.append(degree)
            center_list.append(np.flip(center))
            slope_list.append(slope)
            rms_list.append(rms)

            # Create the center_str for file naming
            center_str = f"{center[1]}_{center[0]}"

            coeff_filename = os.path.join(data_directory,
                                          f"{fit_coeff_path}center_{center_str}_{degree}deg_{n_frames}frames_noframe1.fits")
            residuals_filename = os.path.join(data_directory,
                                              f"{residuals_base_path}center_{center_str}_{degree}deg_{n_frames}frames_noframe1.fits")

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

    # Create the results DataFrame
    results_df = pd.DataFrame({
        'Degree': degree_list,
        'Center': center_list,
        'Slope': slope_list,
        'RMS': rms_list
    })

    # Save the results to a CSV file
    csv_filename = os.path.join(data_directory, 'results_vs_degree_trial.csv')
    results_df.to_csv(csv_filename, index=True)
    print(f"Results saved to {csv_filename}")

    print(results_df)


def plot_superpixels_centers(image_size, centers, size=256):
    fig, ax = plt.subplots()

    # Create an empty image
    ax.imshow(np.zeros(image_size), cmap='gray', origin='lower')

    # Plot centers as dots and rectangles
    half_size = size // 2
    for center in centers:
        # Plot center as a red dot
        ax.plot(center[1], center[0], 'ro', markersize=10)  # Red dots for centers

        # Annotate the center with its coordinates
        ax.text(center[1] + half_size + 10, center[0] + half_size + 10,
                f'({center[1]}, {center[0]})', color='white', fontsize=12,
                verticalalignment='bottom', horizontalalignment='left')

        # Define the rectangle
        rect = patches.Rectangle(
            (center[1] - half_size, center[0] - half_size),
            size, size,
            linewidth=1, edgecolor='r', facecolor='none'
        )
        # Add the rectangle to the plot
        ax.add_patch(rect)

    ax.set_title('Superpixel Centers')
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.grid(True)
    plt.show()


# Parameters
centers = [(500, 2048), (2048, 2048), (3500, 2048)]  # Centers of the superpixels
degrees = np.linspace(1, 10, 10)
generate_fit_cube(degrees, centers)
