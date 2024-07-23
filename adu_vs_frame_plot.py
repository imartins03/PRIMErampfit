from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

y_cube_path = r'D:\NLC\C1\y_cube_100im.fits'
n_frames = 100
super_pixel_size = (256, 256)

#d the y_cube data Loa
y_cube = fits.getdata(y_cube_path)[:n_frames]

#dimensions of the cube
_, y_dim, x_dim = y_cube.shape

y_center, x_center = y_dim//2, x_dim//2
y_start, y_end = y_center-super_pixel_size[0]//2, y_center+super_pixel_size[0]//2   #super_pixel_size[0] is the height so this is taking the center - half the height as the start
x_start, x_end = x_center - super_pixel_size[1] // 2, x_center + super_pixel_size[1] // 2 #super_pixel_size[1] is the width

#extract super pixel for all frames
super_pixel_cube = y_cube[:, y_start:y_end, x_start:x_end]

# Calculate the median ADU for each frame
median_values = [np.median(super_pixel_cube[frame]) for frame in range(n_frames)]

# Plot the median ADU values
plt.plot(range(n_frames), median_values, marker='o')
plt.xlabel('Frame Number')
plt.ylabel('Median ADU Value')
plt.title('Median ADU Value as a Function of Frame Number')
plt.grid(True)
plt.show()
