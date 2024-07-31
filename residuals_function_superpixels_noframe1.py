from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd

# Paths
super_bias_path = 'IRRC_calfiles\\super_biasC1.fits.ramp.20231012'
calFile = r'IRRC_calfiles\irrc_weights_C1.h5'
maskFile_path = r'IRRC_calfiles\C1_bad_ref_pix_mask.fits'
y_cube_path = r'D:\NLC\C1\y_cube_100.fits'

fit_cube_path = r'D:\NLC\C1\dif_degrees_test\fit_cube_poly_8deg.fits'
fit_coeff_path = r'D:\NLC\C1\dif_degrees_test\fit_coeff_poly_8deg.fits'
residuals_cube_path = r'D:\NLC\C1\dif_degrees_test\residuals_poly_8deg.fits'
stat_table = r'D:\NLC\C1\dif_degrees_test\frame_statistics_poly_8deg.csv'

n_frames = 100
superpixel_centers = [(2048, 2048), (3072, 2048), (1024, 2048)]
superpixel_size = 256


def calculate_residuals():
    y_cube = fits.getdata(y_cube_path)[1:n_frames]  # Load y_cube data
    y_cube = y_cube[:, 0, :, :]  # Take out the second dimension

    fit_cube = fits.getdata(fit_cube_path)[1:n_frames]  # Load fit cube data
    residuals_cube = y_cube - fit_cube  # Calculate residuals

    # Save the residuals cube
    fits.writeto(residuals_cube_path, residuals_cube, overwrite=True)

    return residuals_cube


def extract_superpixel_data(residuals_cube, center, size):
    half_size = size // 2
    y_start = center[0] - half_size
    y_end = center[0] + half_size
    x_start = center[1] - half_size
    x_end = center[1] + half_size

    # Extract the superpixel region for all frames
    superpixel_data = residuals_cube[:, y_start:y_end, x_start:x_end]
    return superpixel_data


def calculate_statistics_for_superpixels(superpixel_centers, size, residuals_cube_path, stat_table):
    residuals_cube = fits.getdata(residuals_cube_path)

    # Prepare lists to store statistics
    frame_num = []
    means = []
    rms_vals = []
    median_vals = []
    std_vals = []
    error = []

    # Iterate over each superpixel
    for center in superpixel_centers:
        superpixel_data = extract_superpixel_data(residuals_cube, center, size)

        # Calculate statistics for each frame in the superpixel data
        for i in range(superpixel_data.shape[0]):
            data = superpixel_data[i]
            means.append(np.mean(data))  # Calculate mean
            rms_vals.append(np.sqrt(np.mean(data ** 2)))  # Calculate RMS of residuals
            median_vals.append(np.median(data))  # Calculate median
            std_vals.append(np.std(data))  # Calculate std of residuals
            error.append(np.mean(np.sqrt(np.mean(data ** 2))) / np.sqrt(len(data)))
            frame_num.append(i)  # Frame number

    # Create a DataFrame and save it
    table = pd.DataFrame({
        'Frame': frame_num,
        'Mean': means,
        'RMS': rms_vals,
        'Median': median_vals,
        'StdDev': std_vals,
        'Error': error
    })
    table.to_csv(stat_table, index=False)  # Save statistics to CSV

    # Additional calculations (Total RMS divided by divisor)
    total_rms_sum = np.sum(table['RMS'])
    length_of_data = len(table)
    divisor = np.sqrt(length_of_data - 1)
    table['TotalRMS'] = total_rms_sum
    table['TotalRMSDivided'] = table['TotalRMS'] / divisor

    # Save updated DataFrame
    table.to_csv(stat_table, index=False)

    return table


def plot_superpixels_centers(image_size, centers, size=256):
    fig, ax = plt.subplots()
    ax.imshow(np.zeros(image_size), cmap='gray', origin='lower')

    half_size = size // 2

    # Plot centers as dots and squares
    for center in centers:
        # Plot the center point
        ax.plot(center[1], center[0], 'ro', markersize=10)  # Red dots for centers

        # Define the square's position
        rect = patches.Rectangle(
            (center[1] - half_size, center[0] - half_size),  # Bottom-left corner
            size,  # Width
            size,  # Height
            linewidth=1,  # Line width
            edgecolor='r',  # Edge color
            facecolor='none'  # No fill color
        )
        # Add the rectangle to the plot
        ax.add_patch(rect)

    ax.set_title('Superpixel Centers')
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.grid(True)
    plt.show()


# Run the calculations and plotting
residuals_cube = calculate_residuals()
statistics_df = calculate_statistics_for_superpixels(superpixel_centers, superpixel_size, residuals_cube_path,
                                                     stat_table)

# Plot the superpixels
plot_superpixels_centers(image_size=(4088, 4088), centers=superpixel_centers, size=superpixel_size)

# Print the DataFrame to verify
print(statistics_df.head())
