import os
from astropy.io import fits
import numpy as np
from irrc.apply_correction import _load_coeffs_from_file, apply_in_memory
from irrc.util import destripe
import matplotlib.pyplot as plt

# Define paths and constants
file_format = r'D:\NLC\C1\{0:08d}C1.fits.fz'
r = (1124972, 1124972 + 4)
file_list = [file_format.format(n) for n in range(*r)]
full_file_list = file_list
super_bias_path = 'IRRC_calfiles\\super_biasC1.fits.ramp.20231012'
calFile = r'IRRC_calfiles\irrc_weights_C1.h5'
maskFile_path = r'IRRC_calfiles\C1_bad_ref_pix_mask.fits'


y_cube_path = 'y_cube.fits'
fit_cube_path = 'fit_cube.fits'
residuals_cube_path = 'residuals.fits'

# Load necessary data
super_bias = fits.getdata(super_bias_path)
maskFile = fits.getdata(maskFile_path)
mask = maskFile > 0
supercpy = super_bias.copy()
supercpy[mask[:, :4096]] = 0

def irrc_correct_frame(dataIn, superbias, calFile, multiThread=True, externalPixelFlags=None, superbias_corrected=True):
    """Returns image cube corrected using IRRC in memory, given calibration file."""
    alpha, gamma, zeta = _load_coeffs_from_file(calFile)
    dataIn_tmp = destripe(dataIn[None, :, :])
    dataIn_tmp[:, :, :4096] -= superbias[None, :, :]
    dataIRRC = apply_in_memory(dataIn_tmp, alpha, gamma, zeta, multiThread, externalPixelFlags, superbias_corrected) + superbias
    return dataIRRC

def irrc_correct_frame_file(filename, superbias, calFile, multiThread=True, externalPixelFlags=None, superbias_corrected=True):
    dataIn = fits.getdata(filename)
    dataIn = dataIn[:, 6:]
    dcpy = dataIn.copy()
    dcpy[mask] = 0
    return irrc_correct_frame(dcpy, superbias, calFile, multiThread, externalPixelFlags, superbias_corrected)

def generate_y_cube(frame_list, superbias, calFile, slc=((4, 4092), (4, 4092))):
    y = []

    for frame in frame_list:
        out_img = irrc_correct_frame_file(frame, superbias, calFile)
        if slc is not None:
            out_img = out_img[:, slc[0][0]:slc[0][1], slc[1][0]:slc[1][1]]
        y.append(out_img)
    y = np.asarray(y)

    fits.writeto(y_cube_path, y, overwrite=True)
    return y

def evaluate_poly_array(coeffs, x_array):
    output_arrays = []
    for x in x_array:
        output_array = np.zeros(coeffs.shape[1])
        for n, coeff in enumerate(coeffs):
            output_array += coeff * (x ** n)
        output_arrays.append(output_array)
    return np.asarray(output_arrays)

def generate_fit_cube(frame_list, degrees=1, saturation=50000):
    y_cube = fits.getdata(y_cube_path)
    x = np.arange(len(frame_list))
    y = y_cube.reshape(x.shape[0], -1)
    coefficients, _ = np.polynomial.legendre.legfit(x, y, degrees)
    fit_coeff = coefficients.reshape(degrees + 1, -1)
    fit_cube = evaluate_poly_array(np.flip(fit_coeff, axis=0), x).reshape(len(frame_list), *y_cube.shape[1:])
    fits.writeto(fit_cube_path, fit_cube, overwrite=True)
    return fit_cube

def calculate_residuals():
    y_cube = fits.getdata(y_cube_path)
    fit_cube = fits.getdata(fit_cube_path)
    residuals_cube = y_cube - fit_cube
    fits.writeto(residuals_cube_path, residuals_cube, overwrite=True)
    return residuals_cube

# Create y_cube if it doesn't exist
if not os.path.exists(y_cube_path):
    y_cube = generate_y_cube(full_file_list, supercpy, calFile)
else:
    y_cube = fits.getdata(y_cube_path)

# Create fit_cube if it doesn't exist
if not os.path.exists(fit_cube_path):
    fit_cube = generate_fit_cube(full_file_list, degrees=1)
else:
    fit_cube = fits.getdata(fit_cube_path)

# Calculate residuals
res = calculate_residuals()

# Plot residuals
std = np.nanstd(res)
for i, frame in enumerate(file_list):
    plt.figure()
    residuals_frame = res[i]
    bins = np.arange(-2 * std, 2 * std, std / 20)
    hist = np.histogram(residuals_frame[np.isfinite(residuals_frame)], bins=bins)
    plt.bar(hist[1][:-1], hist[0], color='blue')
    plt.title(f'Histogram of Residuals for Frame {i + 1}')
    plt.xlabel('Residual Value')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()
