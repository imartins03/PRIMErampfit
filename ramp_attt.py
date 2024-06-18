import os
import glob
from astropy.io import fits
import numpy as np
from irrc.apply_correction import _load_coeffs_from_file, apply_in_memory
from irrc.util import destripe
import matplotlib.pyplot as plt
import pandas as pd
from numpy.polynomial.legendre import Legendre

# File format and range for frame files
file_format = r'D:\NLC\C1\{0:08d}C1.fits.fz'
r = (1124972, 1124972+4)

file_list = [file_format.format(n) for n in range(*r)]
full_file_list = file_list

# Load super bias and mask files
super_bias = fits.getdata('IRRC_calfiles\\super_biasC1.fits.ramp.20231012')
calFile = r'IRRC_calfiles\irrc_weights_C1.h5'
maskFile = fits.getdata(r'IRRC_calfiles\C1_bad_ref_pix_mask.fits')
mask = maskFile > 0

supercpy = super_bias.copy()
supercpy[mask[:, :4096]] = 0

def irrc_correct_frame(dataIn, superbias, calFile, multiThread=True, externalPixelFlags=None, superbias_corrected=True):
    alpha, gamma, zeta = _load_coeffs_from_file(calFile)
    dataIn_tmp = destripe(dataIn[None, :, :])
    dataIn_tmp[:, :, :4096] = dataIn_tmp[:, :, :4096] - superbias[None, :, :]
    dataIRRC = apply_in_memory(dataIn_tmp, alpha, gamma, zeta, multiThread, externalPixelFlags, superbias_corrected) + superbias
    return dataIRRC

def irrc_correct_frame_file(filename, superbias, calFile, multiThread=True, externalPixelFlags=None, superbias_corrected=True):
    corrected_file = filename.replace('.fits.fz', '.corrected.fits')
    # if os.path.exists(corrected_file):
    #     return
    existing_file = fits.getdata(corrected_file)
    # else:
    # dataIn = fits.getdata(filename)
    # dataIn = dataIn[:, 6:]
    # dcpy = dataIn.copy()
    # dcpy[mask] = 0
    # corrected_data = irrc_correct_frame(dcpy, superbias, calFile, multiThread, externalPixelFlags, superbias_corrected)
    # #     fits.writeto(corrected_file, corrected_data, overwrite=True)
    return existing_file
    #corrected_data)

# def irrc_correct_frame_file(filename, superbias, calFile, multiThread=True, externalPixelFlags=None, superbias_corrected=True):
#     corrected_file = filename.replace('.fits.fz', '.corrected.fits')
#     dataIn = fits.getdata(filename)
#     dataIn = dataIn[:, 6:]
#     dcpy = dataIn.copy()
#     dcpy[mask] = 0
#     corrected_data = irrc_correct_frame(dcpy, superbias, calFile, multiThread, externalPixelFlags, superbias_corrected)
#     fits.writeto(corrected_file, corrected_data, overwrite=True)
#     return corrected_data

def evaluate_poly_array(coeffs, x_array, poly_type='legendre'):
    output_arrays = []
    for x in x_array:
        if poly_type == 'power':
            output_array = np.zeros(coeffs.shape[1])
            for n, coeff in enumerate(coeffs):
                output_array += coeff*(x**n)
        elif poly_type == 'legendre':
            output_array = np.zeros(coeffs.shape[1])
            for j in range(coeffs.shape[1]):
                leg = Legendre(coeffs[:, j])
                output_array[j] = leg(x)
        output_arrays.append(output_array)
    return np.asarray(output_arrays)

def val(coeff_data_file, frame_list, degrees=1):
    coeff_data = fits.getdata(coeff_data_file)
    num_coeffs = coeff_data.shape[0]
    fit_coeff = coeff_data.reshape(num_coeffs, -1)
    val_fit = evaluate_poly_array(np.flip(fit_coeff, axis=0), np.arange(len(frame_list)), poly_type='legendre')
    return val_fit.reshape(len(frame_list), 4088, 4088)

def compute_y_cb(frame_list, superbias, calFile, mask, degrees=1, saturation=50000):
    y_cb = []
    for frame in frame_list:
        corrected_data = irrc_correct_frame_file(frame, superbias, calFile)
        y_cb.append(corrected_data)
    y_cb = np.asarray(y_cb)
    return y_cb

def residuals(y_cb, fit_cube):
    residuals_cube = y_cb - fit_cube
    return residuals_cube

y_cb_file = full_file_list[0].replace('.fits.fz', '.y_cb.fits')
y_cb = fits.getdata(y_cb_file)

coeff_data_file = full_file_list[0].replace('.fits.fz', f'.img_cb_3deg.fits')
fit_cube = val(coeff_data_file, full_file_list, degrees=3)

# Compute residuals
res = residuals(y_cb, fit_cube)
std = np.nanstd(res)
fits.writeto(r'D:\NLC\C1\residuals_cb.fits', res, overwrite=True)

# Plot histograms for residuals of each frame
for i, frame in enumerate(file_list):
    plt.figure()
    residuals_frame = res[i]
    bins = np.arange(-2*std, 2*std, std/20)
    hist = np.histogram(residuals_frame[np.isfinite(residuals_frame)], bins=bins)
    plt.bar(hist[1][:-1], hist[0], color='blue')
    plt.title(f'Histogram of Residuals for Frame {i+1}')
    plt.xlabel('Residual Value')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()

# Calculate and save means and RMS values for each frame
means = []
rms_vals = []
for frame in file_list:
    data = fits.getdata(frame)
    means.append(np.mean(data))
    rms_vals.append(np.sqrt(np.mean(data**2)))

# Save means and RMS values to CSV
table = pd.DataFrame({'Mean': means, 'RMS': rms_vals})
table.to_csv(r'D:\NLC\C1\frame_statistics.csv', index=False)

# res.shape[0] gives the number of elements along first dimension, corresponds number of frames
for i, frame in enumerate(file_list):
    plt.figure()
    residuals_frame = res[i]
    bins = np.arange(-2*std, 2*std, std/20)
    hist = np.histogram(residuals_frame[np.isfinite(residuals_frame)], bins=bins)
    plt.bar(hist[1][:-1], hist[0], color='blue')
    plt.title(f'Histogram of Residuals for Frame {i+1}')
    plt.xlabel('Residual Value')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()
