from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os


n_frames=239
# Directory where the data files are stored
data_directory = f'F:\laguerrefit\{n_frames}_frames_unweighted'

# Initialize a list to store median residual values for each degree
median_residuals_all_deg = []

# Loop over degrees from 1 to 10
for degree in range(1, 11):

    # Construct file paths
    residuals_cube_path = os.path.join(data_directory, f'res_cube_lag_{n_frames}frames_{degree}deg_final_vst_unweighted.fits')

    # Load the residuals data
    residuals_cube = fits.getdata(residuals_cube_path)


    # List to store median values for the current degree
    median_residuals = []

    # Iterate over each frame
    for frame_n in range(residuals_cube.shape[0]):
        residuals_frame = residuals_cube[frame_n]

        # Compute the median value of the residuals frame
        median_value = np.median(residuals_frame)
        median_residuals.append(median_value)

    # Plot the median residuals as a function of frame number
    plt.figure()
    plt.plot(np.linspace(1,residuals_cube.shape[0],residuals_cube.shape[0]), median_residuals, marker='o', linestyle='-', color='black')
    plt.title(f'Median Residuals ({degree} Degree) as a Function of Frame Number')
    plt.xlabel('Frame Number')
    plt.ylabel('Median Residual Value')


    plt.grid(True)

    # Save the plot
    plot_filename = os.path.join(data_directory, f'median_residuals_plot_lag_{n_frames}frames_{degree}deg_unweighted.png')
    plt.savefig(plot_filename)
    # plt.show()
    plt.close()  # Close the figure to release memory

    # Append median residual values for this degree to the list
    median_residuals_all_deg.append(median_residuals)

print("Analysis complete. Plots saved for each degree.")
