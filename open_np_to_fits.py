# import os
# import numpy as np
# from astropy.io import fits
#
# # Directory containing the .npy files
# np_file_dir = r'F:\legfit'
#
# # Get list of .npy files in the directory
# np_files = [os.path.join(np_file_dir, f) for f in os.listdir(np_file_dir) if f.endswith('.npy')]
#
# def save_coeff_im(degree):
#     for f in np_files:
#         print(f"Processing file: {f}")
#         arr = np.load(f)
#
#         # Print the original shape and size of the array
#         print(f"Original array shape: {arr.shape}")
#         print(f"Original array size: {arr.size}")
#
#         # Check if the array can be reshaped to (degree + 1, 1022, 4088)
#         expected_shape = (arr.shape[0], 1022, 4088)
#
#         arr = arr.reshape(expected_shape)
#
#
#         # Save as .fits file
#         fits_file = f'fit_cube_leg_239frames_{degree}deg_final_vst_unweighted'
#         fits.writeto(fits_file, arr, overwrite=True)
#         print(f"Saved .fits file: {fits_file}")
#
#
# # Process files for degrees from 1 to 10
# for degree in range(1, 11):
#     save_coeff_im(degree)

import os
import numpy as np
from astropy.io import fits

# Directory containing the .npy files
np_file_dir = r'F:\legfit'

# Get list of .npy files in the directory
np_files = [os.path.join(np_file_dir, f) for f in os.listdir(np_file_dir) if f.startswith('coefficients')]

def save_coeff_im(degree):
    for f in np_files:
        print(f"Processing file: {f}")
        arr = np.load(f)

        # Print the original shape and size of the array
        print(f"Original array shape: {arr.shape}")
        print(f"Original array size: {arr.size}")

        # Check if the array can be reshaped to (degree + 1, 1022, 4088)
        expected_shape = (arr.shape[0], 1022, 4088)

        arr = arr.reshape(expected_shape)


        # Save as .fits file
        fits_file = f'fit_cube_leg_239frames_{degree}deg_final_vst_unweighted'
        fits.writeto(fits_file, arr, overwrite=True)
        print(f"Saved .fits file: {fits_file}")


# Process files for degrees from 1 to 10
for degree in range(1, 11):
    save_coeff_im(degree)








