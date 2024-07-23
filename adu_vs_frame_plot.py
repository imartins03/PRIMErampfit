# from astropy.io import fits
# import numpy as np
# import matplotlib.pyplot as plt
#
# y_cube_path = r'D:\NLC\C1\y_cube_100im.fits'
# n_frames = 100
# super_pixel_size = (256, 256)
#
# #d the y_cube data Loa
# y_cube = fits.getdata(y_cube_path)[:n_frames]
#
# #dimensions of the cube
# _, y_dim, x_dim = y_cube.shape
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



import numpy as np
import matplotlib.pyplot as plt

# Assuming 'cube' is a 4D numpy array of shape (100, 4088, 4088)
# with dimensions [frames, height, width]

def extract_super_pixel(image, size=256):
    center = image.shape[0] // 2
    half_size = size // 2
    return image[center - half_size:center + half_size, center - half_size:center + half_size]

def find_median_adus(cube, super_pixel_size=256):
    median_values = []
    for image in cube:
        super_pixel = extract_super_pixel(image, super_pixel_size)
        median_value = np.median(super_pixel)
        median_values.append(median_value)
    return median_values

# Load your 4D cube here
# cube = ...

median_values = find_median_adus(cube)

# Plot median ADU values
plt.plot(range(1, len(median_values) + 1), median_values, marker='o')
plt.xlabel('Frame Number')
plt.ylabel('Median ADU Value')
plt.title('Median ADU Value vs. Frame Number')
plt.show()
  #hey