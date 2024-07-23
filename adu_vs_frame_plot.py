# from astropy.io import fits
# import numpy as np
# import matplotlib.pyplot as plt
#
# y_cube_path = r'D:\NLC\C1\y_cube_4.fits'
# n_frames = 4
# super_pixel_size = (256, 256)
#
# #d the y_cube data Loa
# y_cube = fits.getdata(y_cube_path)[:n_frames]
#
# #dimensions of the cube
# y = y_cube.reshape(y_cube.shape[0],y_cube.shape[2], y_cube.shape[3])
# _, y_dim, x_dim = y.shape
#
#
# y_center, x_center = y_dim//2, x_dim//2
# y_start, y_end = y_center-super_pixel_size[0]//2, y_center+super_pixel_size[0]//2   #super_pixel_size[0] is the height so this is taking the center - half the height as the start
# x_start, x_end = x_center - super_pixel_size[1] // 2, x_center + super_pixel_size[1] // 2 #super_pixel_size[1] is the width
#
# #extract super pixel for all frames
# super_pixel_cube = y_cube[:, y_start:y_end, x_start:x_end]
#
# # Calculate the median ADU for each frame
# median_values = [np.median(super_pixel_cube[frame]) for frame in range(n_frames)]
#
# # Plot the median ADU values
# plt.plot(range(n_frames), median_values, marker='o')
# plt.xlabel('Frame Number')
# plt.ylabel('Median ADU Value')
# plt.title('Median ADU Value as a Function of Frame Number')
# plt.grid(True)
# plt.show()

#keep this stuff
# import numpy as np
# import matplotlib.pyplot as plt

# Assuming 'cube' is a 4D numpy array of shape (100, 4088, 4088)
# with dimensions [frames, height, width]

# def extract_super_pixel(image, size=256):
#     center = image.shape[0] // 2
#     half_size = size // 2
#     return image[center - half_size:center + half_size, center - half_size:center + half_size]
#
# def find_median_adus(cube, super_pixel_size=256):
#     median_values = []
#     for image in cube:
#         super_pixel = extract_super_pixel(image, super_pixel_size)
#         median_value = np.median(super_pixel)
#         median_values.append(median_value)
#     return median_values
#
# # Load your 4D cube here
# # cube = ...
#
# median_values = find_median_adus(cube)
#
# # Plot median ADU values
# plt.plot(range(1, len(median_values) + 1), median_values, marker='o')
# plt.xlabel('Frame Number')
# plt.ylabel('Median ADU Value')
# plt.title('Median ADU Value vs. Frame Number')
# plt.show()
#   #hey
#
# from astropy.io import fits
# import numpy as np
# import matplotlib.pyplot as plt
#
# # Path to y_cube FITS file
# y_cube_path = r'D:\NLC\C1\y_cube_100im.fits'
#
# # Load y_cube data
# y_cube = fits.getdata(y_cube_path)
# y = y_cube.reshape(y_cube.shape[0], y_cube.shape[2], y_cube.shape[3])
#
# # Define the center pixel location (assuming the center is at the middle of the frame)
# center_x, center_y = y.shape[1] // 2, y.shape[2] // 2
#
# # Initialize lists to store center pixel values and frame numbers
# center_pixel_values = []
# frame_numbers = np.arange(y_cube.shape[0])
#
# # Iterate over each frame and extract the center pixel value
# for frame_index in range(y_cube.shape[0]):
#     frame = y_cube[frame_index, 0]  # Select the first element along the second axis
#     # Extract the value of the center pixel
#     center_pixel_value = frame[center_x +30, center_y +50]
#     center_pixel_values.append(center_pixel_value)
#
# # Plot the center pixel values as a function of frame number
# plt.figure(figsize=(10, 6))
# plt.plot(frame_numbers, center_pixel_values, marker='o', linestyle='-', color='b')
# plt.title('Center Pixel Value as a Function of Frame Number')
# plt.xlabel('Frame Number')
# plt.ylabel('Center Pixel Value')
# plt.grid(True)
# plt.savefig(r'D:\NLC\C1\center_pixel_values_plot.png')  # Save the plot
# plt.ylim(0,50000)
# plt.show()

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

residuals_cube_path = r'D:\NLC\C1\residuals.fits'
y = fits.getdata(residuals_cube_path)
center_x, center_y = y.shape[1]//2, y.shape[2]//2
superpix_size = 256
half_size = superpix_size//2

# Initialize lists to store median values and frame numbers
median_values = []
frame_numbers = np.arange(y.shape[0])

# Iterate over each frame
for frame_n in range(y.shape[0]):
    frame = y[frame_n, 0]

    # Extract 256x256 square from the center
    superpix = frame[center_x - half_size:center_x + half_size,
                     center_y - half_size:center_y + half_size]

    median_value = np.median(superpix)
    median_values.append(median_value)

# Plot the median values as a function of frame number
plt.plot(frame_numbers, median_values, marker='o', linestyle='-', color='black')
plt.title('Median of 256x256 super pixel from Center as a Function of Frame Number')
plt.xlabel('Frame Number')
plt.ylabel('Median Value of 256x256 superpixel')
plt.grid(True)
plt.savefig(r'D:\NLC\C1\median_superpix_values_plot.png')  # Save the plot
plt.show()

#save shisdids    bkjbjbk
