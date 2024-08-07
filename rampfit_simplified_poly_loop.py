from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#w
# Definition of paths
super_bias_path = 'IRRC_calfiles\\super_biasC1.fits.ramp.20231012'
calFile = r'IRRC_calfiles\irrc_weights_C1.h5'
maskFile_path = r'IRRC_calfiles\C1_bad_ref_pix_mask.fits'
y_cube_path = r'D:\NLC\C1\y_cube_500.fits'
fit_cube_path_template = r'F:\leftover_C1_dif_degrees_test_rampfit\239_frames\fit_cube_poly_{degree}deg_239frames_noframe1.fits'
fit_coeff_path_template = r'F:\leftover_C1_dif_degrees_test_rampfit\239_frames\fit_coeff_poly_{degree}deg_239frames_noframe1.fits'
residuals_cube_path_template = r'F:\leftover_C1_dif_degrees_test_rampfit\239_frames\residuals_poly_{degree}deg_239frames_noframe1.fits'

# Getting data from existing files
super_bias = fits.getdata(super_bias_path)  # Load super bias data
maskFile = fits.getdata(maskFile_path)  # Load mask file data
mask = maskFile > 0  # Create mask for bad pixels
supercpy = super_bias.copy()  # Create a copy of the super bias
supercpy[mask[:, :4096]] = 0  # Apply mask to the super bias copy

def evaluate_poly_array(coeffs, a_array, poly_type='power'):
    # Function to evaluate polynomial arrays
    output_arrays = []
    for a in a_array:  # Loop over input values
        if poly_type == 'power':  # Only 'power' type is considered here
            output_array = np.zeros(coeffs.shape[1])  # Initialize output array
            for n, coeff in enumerate(coeffs):  # Loop over coefficients
                output_array += coeff * (a ** n)  # Calculate polynomial value
            output_arrays.append(output_array)  # Append result to list
    return np.asarray(output_arrays)  # Convert list to numpy array

def generate_fit_cube(degree, saturation=50000, n_frames=239):
    y_cube = fits.getdata(y_cube_path)  # Load y_cube data
    x = y_cube.shape[0]  # x is the dimension of the data cube (number of frames)

    y = y_cube.reshape(x, 4088, 4088)
    if n_frames is not None:
        y = y[1:n_frames, :, :]  # Use specified number of frames
        x = y.shape[0]

    y_cube = y
    y = y_cube.reshape(x, -1)  # Reshape y_cube for fitting

    time = np.arange(len(y), dtype=np.double)

    # Generate array for fitting, time in units of frames
    coefficients, _ = np.polyfit(time, y, degree, cov=True)  # Fit polynomial

    # Reshape coefficients and save
    fit_coeff = coefficients.reshape(degree + 1, 4088, 4088)


    fit_coeff_path = fit_coeff_path_template.format(degree=degree)
    fits.writeto(fit_coeff_path, fit_coeff, overwrite=True)

    fit_cube = evaluate_poly_array(np.flip(coefficients, axis=0), time)  # Evaluate polynomial array

    fit_cube = fit_cube.reshape(len(time), 4088, 4088)  # Reshape fit cube
    fit_cube_path = fit_cube_path_template.format(degree=degree)
    fits.writeto(fit_cube_path, fit_cube, overwrite=True)

    return fit_cube

# Loop through polynomial degrees from 1 to 10
saturation = 50000  # Currently not used
for degree in range(1, 11):
    print(f"Processing degree {degree}")
    generate_fit_cube(degree, saturation)
