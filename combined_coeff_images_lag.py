from astropy.io import fits
import numpy as np
import os

# Directory containing the .fits files
input_dir = r'F:\laguerrefit\239_frames_unweighted'
output_dir = r'F:\laguerrefit\239_frames_unweighted\coeff_im_final'


def reshape_and_save_images(degree):
    row_images = []

    # Loop through each row (0 to 3)
    for row in range(4):
        # Construct the file name for the given degree and row
        file_name = f'coefficients_239frames_degree{degree}_row{row}_weighted.fits'
        file_path = os.path.join(input_dir, file_name)

        # Load the FITS file
        with fits.open(file_path) as hdul:
            arr = hdul[0].data

        # Print the original shape and size of the array
        print(f"Original array shape: {arr.shape}")
        print(f"Original array size: {arr.size}")

        # Reshape the array to (degree + 1, 1022, 4088) if needed
        reshaped_arr = arr.reshape(degree + 1, 1022, 4088)
        row_images.append(reshaped_arr)
#
        # Save the reshaped image back to disk
        reshaped_file_path = os.path.join(input_dir, f'coefficients_lag_239frames_{degree}deg_final_vst_weighted.fits')
        fits.writeto(reshaped_file_path, reshaped_arr, overwrite=True)
        print(f"Saved reshaped .fits file: {reshaped_file_path}")

    # Stack along the row axis (axis=1), which combines the row images
    combined_image = np.concatenate(row_images, axis=1)

    # Check if the combined shape is as expected
    expected_shape = (degree + 1, 4088, 4088)
    assert combined_image.shape == expected_shape, f"Unexpected shape: {combined_image.shape}, expected: {expected_shape}"

    # Save the final combined image
    final_file_path = os.path.join(output_dir, f'coefficients_lag_239frames_{degree}deg_final_vst_weighted.fits')
    fits.writeto(final_file_path, combined_image, overwrite=True)
    print(f"Saved final .fits file: {final_file_path}")

