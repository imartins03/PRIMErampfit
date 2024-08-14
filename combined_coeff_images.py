from astropy.io import fits
import numpy as np
import os

# Directory containing the .npy files
output_dir = r'F:\legfit\239_frames_unweighted'

np_files = [os.path.join(output_dir, f) for f in os.listdir(output_dir) if f.startswith('coefficients')]

def save_full_coeff_im(degree):
    for f in np_files:
        print(f"Processing file: {f}")
        arr = np.load(f)

        # Print the original shape and size of the array
        print(f"Original array shape: {arr.shape}")
        print(f"Original array size: {arr.size}")

        # # Check if the array can be reshaped to (degree + 1, 1022, 4088)
        # expected_shape = (arr.shape[0], 1022, 4088)
        #
        # arr = arr.reshape(expected_shape)

        # Save final stacked fit cube
        fits_file_path = os.path.join(output_dir, f'coefficients_leg_239frames_{degree}deg_final_vst_unweighted')

        fits.writeto(fits_file_path, arr, overwrite=True)
        print(f"Saved .fits file: {fits_file_path}")


# Process files for degrees from 1 to 10
for degree in range(1, 11):
    save_full_coeff_im(degree)