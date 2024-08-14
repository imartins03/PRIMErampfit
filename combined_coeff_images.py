from astropy.io import fits
import numpy as np
import os

# Directory containing the .fits files
input_dir = r'F:\legfit\239_frames_unweighted'
output_dir = r'F:\legfit\239_frames_unweighted'


def combine_and_save_coeff_images(degree):
    # Prepare a list to hold the image arrays for each row
    row_images = []

    # Loop through each row (0 to 3)
    for row in range(4):
        # Construct the file name for the given degree and row
        file_name = f'Coefficients_239frames_degree{degree}_row{row}_unweighted.fits'
        file_path = os.path.join(input_dir, file_name)

        print(f"Processing file: {file_path}")

        # Load the FITS file
        with fits.open(file_path) as hdul:
            arr = hdul[0].data

        # Print the original shape and size of the array
        print(f"Original array shape: {arr.shape}")
        print(f"Original array size: {arr.size}")

        # Append to the list
        row_images.append(arr)

    # Stack along the row axis (axis=1), which is the second dimension
    combined_image = np.concatenate(row_images, axis=1)

    # Check if the combined shape is as expected
    expected_shape = (degree + 1, 4088, 4088)
    assert combined_image.shape == expected_shape, f"Unexpected shape: {combined_image.shape}, expected: {expected_shape}"

    # Save the final combined image
    fits_file_path = os.path.join(output_dir, f'coefficients_leg_239frames_{degree}deg_final_vst_unweighted.fits')
    fits.writeto(fits_file_path, combined_image, overwrite=True)
    print(f"Saved .fits file: {fits_file_path}")


# Process files for degrees from 1 to 10
for degree in range(1, 11):
    combine_and_save_coeff_images(degree)
