from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

# Paths
super_bias_path = 'IRRC_calfiles\\super_biasC1.fits.ramp.20231012'
maskFile_path = r'IRRC_calfiles\C1_bad_ref_pix_mask.fits'
y_cube_path = r'D:\NLC\C1\y_cube_100.fits'
data_directory = r'D:\NLC\C1\superpix'
#
# Load data
super_bias = fits.getdata(super_bias_path)  # Load super bias data
maskFile = fits.getdata(maskFile_path)  # Load mask file data
mask = maskFile > 0  # Create mask for bad pixels
supercpy = super_bias.copy()  # Create a copy of the super bias
supercpy[mask[:, :4096]] = 0  # Apply mask to the super bias copy
n_frames = 100

def generate_fit_cube(degrees, centers):
    y_cube = fits.getdata(y_cube_path)[1:n_frames]
    x, _, y, z = y_cube.shape
    # y = y_cube.reshape(x, 4088, 4088)
    y_cube = y_cube[:,0,:,:]

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

            # plt.figure()

            plt.plot(time, residuals, label=f'Center {np.flip(center)}')
            print(f'Center {center} Degree {degree}')


        plt.title(f'Residuals vs. Frame Number for Degree {degree}')
        plt.xlabel('Frame Number')
        plt.ylabel('Residuals')
        plt.legend()
        plt.grid(True)

        plot_filename = os.path.join(data_directory, f'median_superpix_plot_poly_{degree}deg.png')
        plt.savefig(plot_filename)

        plt.show()
        plt.close()

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
centers = [(500, 2048),(2048, 2048), (3500, 2048)]  # Centers of the superpixels
degrees = np.linspace(1,10,10)
generate_fit_cube(degrees, centers)
