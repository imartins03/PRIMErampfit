from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

n_frames=22

y_cube_path = r"G:\20240820.rimas.0010-0109.HK.mean.fits"

# fit_cube_path_template = r'D:\NLC\C1\dif_degrees_test\100_frames\fit_cube_poly_{degree}deg_{n_frames}frames_noframe1_weights.fits'
# fit_coeff_path_template = r'D:\NLC\C1\dif_degrees_test\100_frames\fit_coeff_poly_{degree}deg_{n_frames}frames_noframe1_weights.fits'
# residuals_cube_path_template = r'F:\leftover_C1_dif_degrees_test_rampfit\239_frames\residuals_poly_{degree}deg_239frames_noframe1.fits'

#239 frame paths
fit_cube_path_template = r"G:\fit_cube_poly_{degree}deg_{n_frames}frames.fits"
fit_coeff_path_template = r"G:\fit_coeff_poly_{degree}deg_{n_frames}frames.fits"
residuals_cube_path_template = r"G:\residuals_poly_{degree}deg_{n_frames}frames.fits"

# def weight_func(t):
#     return (50**2)/(50**2 + 306*t)

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

def generate_fit_cube(degree, n_frames=22):
    y_cube = fits.getdata(y_cube_path)  # Load y_cube data
    x = y_cube.shape[0]  # x is the dimension of the data cube (number of frames)

    y = y_cube
    if n_frames is not None:
        y = y[:n_frames, :, :]  # Use specified number of frames
        x = y.shape[0]

    y = y_cube.reshape(x, -1)  # Reshape y_cube for fitting

    time = np.arange(len(y), dtype=np.double)

    # weights = weight_func(time+2)
    # Generate array for fitting, time in units of frames
    # coefficients, _ = np.polyfit(time, y, degree, w=weights, cov=True)  # Fit polynomial
    coefficients, _ = np.polyfit(time, y, degree, cov=True)  # Fit polynomial
    # Reshape coefficients and save
    fit_coeff = coefficients.reshape(degree + 1, y_cube.shape[1], y_cube.shape[2])


    fit_coeff_path = fit_coeff_path_template.format(degree=degree,n_frames=n_frames)
    fits.writeto(fit_coeff_path, fit_coeff, overwrite=True)

    fit_cube = evaluate_poly_array(np.flip(coefficients, axis=0), time)  # Evaluate polynomial array

    fit_cube = fit_cube.reshape(len(time), y_cube.shape[1], y_cube.shape[2])  # Reshape fit cube
    fit_cube_path = fit_cube_path_template.format(degree=degree,n_frames=n_frames)
    fits.writeto(fit_cube_path, fit_cube, overwrite=True)

    return fit_cube


for degree in range(1, 11):
    print(f"Processing degree {degree}")
    generate_fit_cube(degree)
