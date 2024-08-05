from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

# Directory where the data files are stored
data_directory = r'F:\leftover_C1_dif_degrees_test_rampfit\239_frames'

# Initialize a list to store slopes and degrees
slopes = []
degrees = list(range(1, 11))

# Loop over all degrees
for degree in degrees:
    # Construct the file path for residuals data
    residuals_cube_path = os.path.join(data_directory, f'residuals_poly_{degree}deg_239frames_noframe1.fits')

    # Load the residuals data
    residuals_cube = fits.getdata(residuals_cube_path)

    median_values = []

    # Iterate over each frame
    for frame_n in range(residuals_cube.shape[0]):
        frame = residuals_cube[frame_n]
        median_value = np.median(frame)
        median_values.append(median_value)

    # Compute the slope of median values
    x_values = np.arange(len(median_values))
    slope, _ = np.polyfit(x_values, median_values, 1)

    # Store the slope for this degree
    slopes.append(slope)

    # Plot the median values as a function of frame number
    # plt.figure()
    # plt.plot(x_values, median_values, marker='o', linestyle='-', color='black')
    # plt.title(f'Median of Whole Image ({degree} deg fit) from Center as a Function of Frame Number (poly)')
    # plt.xlabel('Frame Number (no frame 1)')
    # plt.ylabel('Median Value of Image')
    # plt.grid(True)
    #
    # plt.annotate(f'Slope: {slope:.4f}', xy=(0.95, 0.95), xycoords='axes fraction', fontsize=10,
    #              verticalalignment='top', horizontalalignment='right', color='red')
    #
    # # Save the plot
    # plot_filename = os.path.join(data_directory, f'median_fullim_plot_poly_{degree}deg_noframe1.png')
    # plt.savefig(plot_filename)
    # plt.show()
    # plt.close()  # Close the figure to release memory

# Create a table of slopes vs. degrees
slope_table = pd.DataFrame({
    'Degree': degrees,
    'Slope': slopes
})

# Display the table
print("Slope vs. Degree Table:")
print(slope_table)

# Save table to CSV
csv_filename = os.path.join(data_directory, 'slope_vs_degree.csv')
slope_table.to_csv(csv_filename, index=True)
print(f"Saved slope vs. degree table to '{csv_filename}'.")
