from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Definition of paths
super_bias_path = 'IRRC_calfiles\\super_biasC1.fits.ramp.20231012'
calFile = r'IRRC_calfiles\irrc_weights_C1.h5'
maskFile_path = r'IRRC_calfiles\C1_bad_ref_pix_mask.fits'
y_cube_path = r'D:\NLC\C1\y_cube_500.fits'

def evaluate_legendre_poly(coeffs, x):
    return np.polynomial.legendre.legval(x, coeffs)

def generate_fit_cube(degree, saturation=50000, n_frames=100):
    y_cube = fits.getdata(y_cube_path)  # Load y_cube data
    print(f"Processing degree {degree}")
    x = y_cube.shape[0]  # x is the dimension of the data cube (number of frames)

    y = y_cube.reshape(x, 4088, 4088)
    if n_frames is not None:
        y = y[1:n_frames, :, :]  # Use specified number of frames
        x = n_frames

    print(y.shape)
    y = y_cube.reshape(x, -1)  # Reshape y_cube for fitting

    time = np.linspace(-1, 1, x, dtype=np.double)  # Legendre time

    # Fit Legendre polynomial
    coefficients = np.polynomial.legendre.legfit(time, y, degree)
    print('Fitting coefficients:', coefficients)

    # Save fit coefficients cube
    fit_coeff = coefficients.reshape(degree + 1, 4088, 4088)  # Reshape coefficients
    fit_coeff_path = f'F:\\legfit\\fit_coeff_leg_{degree}deg.fits'
    fits.writeto(fit_coeff_path, fit_coeff, overwrite=True)  # Save coefficients
    print(fit_coeff.shape)

    # Generate and save fit cube
    fit_cube = evaluate_legendre_poly(coefficients, time)
    fit_cube = fit_cube.reshape(len(time), 4088, 4088)  # Reshape fit cube
    fit_cube_path = f'F:\\legfit\\fit_cube_leg_{degree}deg_noframe1.fits'
    fits.writeto(fit_cube_path, fit_cube, overwrite=True)  # Save fit cube
    print('Fit cube shape:', fit_cube.shape)

# Loop through degrees from 1 to 10
for degree in range(1, 11):
    generate_fit_cube(degree)
