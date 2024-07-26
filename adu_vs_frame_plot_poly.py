from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

# Directory where the data files are stored
data_directory = r'D:\NLC\C1\dif_degrees_test'
#
# Initialize lists to store median values and frame numbers
# median_values_all_deg = []
#
# # Loop over all 10 degrees
# for degree in range(1, 11):
#     # Construct the file path for residuals data
#     residuals_cube_path = os.path.join(data_directory, f'residuals_poly_{degree}deg.fits')
#
#     # Load the data
#     y = fits.getdata(residuals_cube_path)
#     center_x, center_y = y.shape[1] // 2, y.shape[2] // 2
#     superpix_size = 256
#     half_size = superpix_size // 2
#
#     median_values = []
#
#     # Iterate over each frame
#     for frame_n in range(y.shape[0]):
#         frame = y[frame_n]
#
#         # Extract 256x256 square from the center
#         superpix = frame[center_x - half_size:center_x + half_size,
#                    center_y - half_size:center_y + half_size]
#
#         median_value = np.median(superpix)
#         median_values.append(median_value)
#
#     # Plot the median values as a function of frame number
#     plt.figure()
#     plt.plot(np.arange(len(median_values)), median_values, marker='o', linestyle='-', color='black')
#     plt.title(f'Median of 256x256 super pixel ({degree} deg fit) from Center as a Function of Frame Number (polyfit)')
#     plt.xlabel('Frame Number')
#     plt.ylabel('Median Value of 256x256 superpixel')
#     plt.grid(True)
#
#     # Save the plot
#     plot_filename = os.path.join(data_directory, f'median_superpix_plot_poly_{degree}deg.png')
#     plt.savefig(plot_filename)
#     plt.show()
#     plt.close()  # Close the figure to release memory
#
#     # Append median values for this degree to the list
#     median_values_all_deg.append(median_values)

#%%
# Initialize lists to store median values and frame numbers
# median_values_all_deg = []
#
# # Loop over all 10 degrees
# for degree in range(1, 11):
#     # Construct the file path for residuals data
#     residuals_cube_path = os.path.join(data_directory, f'residuals_poly_{degree}deg.fits')
#
#     # Load the data
#     y = fits.getdata(residuals_cube_path)
#     center_x, center_y = y.shape[1] // 2, y.shape[2] // 2
#     superpix_size = 256
#     half_size = superpix_size // 2
#
#     median_values_center = []
#     median_values_above = []
#     median_values_below = []
#
#     # Define the vertical distance for the above and below superpixels
#     vertical_distance = (center_y - half_size) // 2
#
#     # Iterate over each frame
#     for frame_n in range(y.shape[0]):
#         frame = y[frame_n]
#
#         # Extract 256x256 square from the center
#         superpix_center = frame[center_x - half_size:center_x + half_size, center_y - half_size:center_y + half_size]
#
#         # Extract 256x256 square above the center
#         superpix_above = frame[center_x - half_size:center_x + half_size,
#                          center_y - half_size - vertical_distance:center_y + half_size - vertical_distance]
#
#         # Extract 256x256 square below the center
#         superpix_below = frame[center_x - half_size:center_x + half_size,
#                          center_y - half_size + vertical_distance:center_y + half_size + vertical_distance]
#
#         median_value_center = np.median(superpix_center)
#         median_value_above = np.median(superpix_above)
#         median_value_below = np.median(superpix_below)
#
#         median_values_center.append(median_value_center)
#         median_values_above.append(median_value_above)
#         median_values_below.append(median_value_below)
#
#     # Plot the median values as a function of frame number
#     plt.figure()
#     plt.plot(np.arange(len(median_values_center)), median_values_center, marker='o', linestyle='-', color='black',
#              label='Center')
#     plt.plot(np.arange(len(median_values_above)), median_values_above, marker='o', linestyle='-', color='red',
#              label='Above Center')
#     plt.plot(np.arange(len(median_values_below)), median_values_below, marker='o', linestyle='-', color='blue',
#              label='Below Center')
#     plt.title(f'Median of 256x256 super pixels ({degree} deg fit) as a Function of Frame Number (polyfit)')
#     plt.xlabel('Frame Number')
#     plt.ylabel('Median Value of 256x256 superpixels')
#     plt.legend()
#     plt.grid(True)
#
#     # Save the plot
#     plot_filename = os.path.join(data_directory, f'median_superpix_plot_poly_{degree}deg.png')
#     plt.savefig(plot_filename)
#     plt.show()
#     plt.close()  # Close the figure to release memory
#
#     # Append median values for this degree to the list
#     median_values_all_deg.append((median_values_center, median_values_above, median_values_below))

#%%

median_values_all_deg = []

# Loop over all 10 degrees
for degree in range(1, 11):
    # Construct the file path for residuals data
    residuals_cube_path = os.path.join(data_directory, f'residuals_poly_{degree}deg.fits')

    # Load the data
    y = fits.getdata(residuals_cube_path)
    # center_x, center_y = y.shape[1] // 2, y.shape[2] // 2
    # superpix_size = 256
    # half_size = superpix_size // 2

    median_values = []

    # Iterate over each frame
    for frame_n in range(y.shape[0]):
        frame = y[frame_n]

        # Extract 256x256 square from the center
        # superpix = frame[center_x - half_size:center_x + half_size,
        #            center_y - half_size:center_y + half_size]

        median_value = np.median(frame)
        median_values.append(median_value)

    # Plot the median values as a function of frame number
    df = pd.read_csv(r'D:\NLC\C1\frame_statistics_poly_10deg.csv')

    err_bar = df['Error'].values

    plt.figure()

    x_values = np.arange(len(median_values))  # Create x-values for frames
    plt.errorbar(x_values, median_values, yerr=err_bar, fmt='o', linestyle='-', color='black')

    if degree <= 2:
        plt.xlim(-.5, 101)
        plt.ylim(auto=True)
    else:
        plt.xlim(-.5, 101)
        plt.ylim(-4.5, 3)

    plt.title(f'Median of whole image ({degree} deg fit) from Center as a Function of Frame Number (poly)')
    plt.xlabel('Frame Number')
    plt.ylabel('Median Value of image')
    plt.grid(True)

    # Save the plot
    plot_filename = os.path.join(data_directory, f'median_fullim_plot_poly_{degree}deg.png')
    plt.savefig(plot_filename)
    plt.show()
    plt.close()  # Close the figure to release memory

    # Append median values for this degree to the list
    median_values_all_deg.append(median_values)