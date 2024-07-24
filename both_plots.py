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

residuals_cube_path = r'D:\NLC\C1\residuals_poly.fits'
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
median_values_poly = []
frame_numbers_poly = np.arange(y.shape[0])

# Iterate over each frame
for frame_n in range(y.shape[0]):
    frame = y[frame_n]

    # Extract 256x256 square from the center
    superpix = frame[center_x - half_size:center_x + half_size,
                     center_y - half_size:center_y + half_size]

    median_value_poly = np.median(superpix)
    median_values_poly.append(median_value_poly)

# Plot the median values as a function of frame number
plt.figure(1)


#save stuff
residuals_cube_path_poly = r'D:\NLC\C1\residuals_leg.fits'
y_leg = fits.getdata(residuals_cube_path_poly)
center_x, center_y = y_leg.shape[1]//2, y_leg.shape[2]//2
superpix_size = 256
half_size = superpix_size//2

# fit_cube_path = r'D:\NLC\C1\fit_cube_poly.fits'
# y = fits.getdata(fit_cube_path)
# center_x, center_y = y.shape[1]//2, y.shape[2]//2# Initialize lists to store median values and frame numbers
# superpix_size = 256median_values = []
# half_size = superpix_size//2frame_numbers = np.arange(y.shape[0])

median_values_leg = []
frame_numbers_leg = np.arange(y_leg.shape[0])

# Iterate over each frame
for frame_n in range(y_leg.shape[0]):
    frame = y_leg[frame_n]

    # Extract 256x256 square from the center
    superpix = frame[center_x - half_size:center_x + half_size,
                     center_y - half_size:center_y + half_size]

    median_value_leg = np.median(superpix)
    median_values_leg.append(median_value_leg)

# Plot the median values as a function of frame number
plt.figure(1)
plt.plot(frame_numbers_leg, median_values_leg, marker='o',linestyle = '-.', color='black',label = 'legendre')
plt.plot(frame_numbers_poly, median_values_poly, linestyle='-', color='red',alpha=0.3,label = 'polynomial')
plt.title('Median of 256x256 super pixel from Center as a Function of Frame Number')
plt.xlabel('Frame Number')
plt.ylabel('Median Value of 256x256 superpixel')
plt.grid(True)
plt.savefig(r'D:\NLC\C1\median_superpix_values_plot_both.png')  # Save the plot
plt.legend()
plt.show()
