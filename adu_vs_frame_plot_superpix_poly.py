from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

data_directory = r'D:\NLC\C1\superpix'

superpixel_centers = [(2048, 2048), (3072, 2048), (500, 2048)]  # Add more centers if needed
median_values_all_deg = []

# Loop over all degrees
for degree in range(1, 11):  # Assuming 10 degrees
    #store median values for each superpixel
    median_values_dict = {center: [] for center in superpixel_centers}

    # Loop over all superpixels
    for center in superpixel_centers:
        # Construct the file path for residuals data
        residuals_cube_path = os.path.join(data_directory, f'residuals_center_{center[1]}_{center[0]}_{degree}deg_noframe1.fits')

        # Load the data
        y = fits.getdata(residuals_cube_path)
        center_x, center_y = y.shape[1] // 2, y.shape[2] // 2
        superpix_size = 256
        half_size = superpix_size // 2

        # Iterate over each frame
        for frame_n in range(y.shape[0]):
            frame = y[frame_n]

            #Extract 256x256 square from the superpixel center, same as previous code
            superpix = frame[center_x - half_size:center_x + half_size, center_y - half_size:center_y + half_size]

            median_value = np.median(superpix)
            median_values_dict[center].append(median_value)

    # Plot the median values as a function of frame number
    plt.figure()
    for center in superpixel_centers:
        plt.plot(np.arange(1, len(median_values_dict[center]) + 1), median_values_dict[center],
                 marker='o', linestyle='-', label=f'Center {center}')

    plt.title(f'Median of 256x256 Superpixels ({degree} deg fit) as a Function of Frame Number')
    plt.xlabel('Frame Number')
    plt.ylabel('Median Value of 256x256 Superpixels')
    plt.legend()
    plt.grid(True)

    # Save the plot
    plot_filename = os.path.join(data_directory, f'median_superpix_plot_poly_{degree}deg.png')
    plt.savefig(plot_filename)
    plt.show()
    plt.close()  # Close the figure to release memory

    # Append median values for this degree to the list
    median_values_all_deg.append(median_values_dict)


