from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# Path to y_cube FITS file
# y_cube_path = r'D:\NLC\C1\y_cube_100im.fits'
# y_cube = fits.getdata(y_cube_path)
# y = y_cube.reshape(y_cube.shape[0], y_cube.shape[2], y_cube.shape[3])
# center_x, center_y = y.shape[1]//2, y.shape[2]//2
# superpix_size = 256
# half_size = superpix_size//2

residuals_cube_path = r'D:\NLC\C1\residuals_leg.fits'
y = fits.getdata(residuals_cube_path)
center_x, center_y = y.shape[1]//2, y.shape[2]//2
superpix_size = 256
half_size = superpix_size//2

# fit_cube_path = r'D:\NLC\C1\fit_cube_leg.fits'
# y = fits.getdata(fit_cube_path)
# center_x, center_y = y.shape[1]//2, y.shape[2]//2
# superpix_size = 256
# half_size = superpix_size//2

# Initialize lists to store median values and frame numbers
median_values = []
frame_numbers = np.arange(y.shape[0])

# Iterate over each frame
for frame_n in range(y.shape[0]):
    frame = y[frame_n]

    # Extract 256x256 square from the center
    superpix = frame[center_x - half_size:center_x + half_size,
                     center_y - half_size:center_y + half_size]

    median_value = np.median(superpix)
    median_values.append(median_value)

# Plot the median values as a function of frame number
plt.figure(1)


#save stuff
residuals_cube_path_poly = r'D:\NLC\C1\residuals_poly.fits'
y2 = fits.getdata(residuals_cube_path_poly)
center_x, center_y = y2.shape[1]//2, y2.shape[2]//2
superpix_size = 256
half_size = superpix_size//2

# fit_cube_path = r'D:\NLC\C1\fit_cube_poly.fits'
# y = fits.getdata(fit_cube_path)
# center_x, center_y = y.shape[1]//2, y.shape[2]//2# Initialize lists to store median values and frame numbers
# superpix_size = 256median_values = []
# half_size = superpix_size//2frame_numbers = np.arange(y.shape[0])

median_values2 = []
frame_numbers2 = np.arange(y2.shape[0])

# Iterate over each frame
for frame_n in range(y2.shape[0]):
    frame = y2[frame_n]

    # Extract 256x256 square from the center
    superpix = frame[center_x - half_size:center_x + half_size,
                     center_y - half_size:center_y + half_size]

    median_value2 = np.median(superpix)
    median_values2.append(median_value2)

# Plot the median values as a function of frame number
plt.figure(1)
plt.plot(frame_numbers2, median_values2, marker='o',linestyle = '-.', color='black')
plt.plot(frame_numbers, median_values, linestyle='-', color='red',alpha=0.3)
plt.title('Median of 256x256 super pixel from Center as a Function of Frame Number')
plt.xlabel('Frame Number')
plt.ylabel('Median Value of 256x256 superpixel')
plt.grid(True)
plt.savefig(r'D:\NLC\C1\median_superpix_values_plot_both.png')  # Save the plot
plt.show()
