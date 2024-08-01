from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

# Directory where the data files are stored
data_directory = r'D:\NLC\C1\dif_degrees_test'

# Initialize a list to store median values across all degrees
median_values_all_deg = []

# Loop over all 10 degrees
for degree in range(1, 11):
    # Construct the file path for residuals data
    residuals_cube_path = os.path.join(data_directory, f'residuals_poly_{degree}deg_noframe1.fits')

    # Load the data
    y = fits.getdata(residuals_cube_path)

    median_values = []

    # Iterate over each frame
    for frame_n in range(y.shape[0]):
        frame = y[frame_n]
        median_value = np.median(frame)
        median_values.append(median_value)

    # Plot the median values as a function of frame number
    plt.figure()
    x_values = np.arange(len(median_values))  # Create x-values for frames
    plt.plot(x_values, median_values, marker='o', linestyle='-', color='black')

    if degree <= 2:
        plt.xlim(-0.5, 101)
        plt.ylim(auto=True)
    else:
        plt.xlim(-0.5, 101)
        plt.ylim(-4.5, 3)

    plt.title(f'Median of Whole Image ({degree} deg fit) from Center as a Function of Frame Number (poly)')
    plt.xlabel('Frame Number (no frame 1)')
    plt.ylabel('Median Value of Image')
    plt.grid(True)

    # Save the plot
    plot_filename = os.path.join(data_directory, f'median_fullim_plot_poly_{degree}deg_noframe1.png')
    plt.savefig(plot_filename)
    plt.show()
    plt.close()  # Close the figure to release memory

    # Append median values for this degree to the list
    median_values_all_deg.append(median_values)

