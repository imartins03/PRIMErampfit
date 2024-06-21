from Demos.BackupRead_BackupWrite import outfile
import math
import sys
import os
import glob
from astropy.io import fits
import numpy as np
from irrc.apply_correction import _load_coeffs_from_file, apply_in_memory
from irrc.util import destripe
import matplotlib.pyplot as plt
import pandas as pd
import scipy
from scipy.optimize import curve_fit
from astropy.modeling.polynomial import Polynomial1D, Legendre1D
from astropy.modeling import fitting
from numpy.polynomial.legendre import Legendre
from numpy.polynomial.polynomial import Polynomial

file_format = r'D:\NLC\C1\{0:08d}C1.fits.fz'
r = (1124972, 1124972+100)

file_list = [file_format.format(n) for n in range(*r)]
full_file_list = file_list

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
    if os.path.exists(corrected_file):
        return fits.getdata(corrected_file)
    else:
        dataIn = fits.getdata(filename)
        dataIn = dataIn[:, 6:]
        dcpy = dataIn.copy()
        dcpy[mask] = 0
        corrected_data = irrc_correct_frame(dcpy, superbias, calFile, multiThread, externalPixelFlags, superbias_corrected)
        fits.writeto(corrected_file, corrected_data, overwrite=True)
        return corrected_data

def get_ramp_slope(frame_list, superbias, calFile, mask, slc=((4, 4092), (4, 4092)), degrees=1, saturation=50000):
    coeff_file = frame_list[0].replace('.fits.fz', f'.img_cb_{degrees}deg.fits')
    if os.path.exists(coeff_file):
        coeff_images = fits.getdata(coeff_file)
    else:
        y = []
        for frame in frame_list:
            out_img = irrc_correct_frame_file(frame, superbias, calFile)
            if slc is not None:
                out_img = out_img[:, slc[0][0]:slc[0][1], slc[1][0]:slc[1][1]]
            saturation_mask = out_img > saturation
            y.append(out_img)
        x = np.arange(len(y))
        y = np.asarray(y)
        y = y.reshape((x.shape[0], -1))
        sat_pix = (y > saturation).astype(float)
        # coefficients, _ = np.polyfit(x, y, degrees, cov=True)
        coefficients, _ = np.polynomial.legendre.legfit(x, y, degrees)
        coefficients[:, sat_pix.sum(axis=0) > 0] = np.nan
        coeff_images = np.stack([coeff.reshape(out_img.shape) for coeff in coefficients], axis=0)
        fits.HDUList([fits.PrimaryHDU(coeff_images)]).writeto(coeff_file, overwrite=True)
    return coeff_images

# def evaluate_poly_array(coeffs, x_array):
#     output_arrays = []
#     for x in x_array:
#         output_array = np.zeros(coeffs.shape[1])
#         for n, coeff in enumerate(coeffs):
#             output_array += coeff * (x ** n)
#         output_arrays.append(output_array)
#     return np.asarray(output_arrays)

#nont poly array i just named it that so it would transfer easily
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

def val(frame_list, superbias, calFile, mask, slc=((4, 4092), (4, 4092)), degrees=1, saturation=50000):
    coeff_data_file = full_file_list[0].replace('.fits.fz', f'.img_cb_{degrees}deg.fits')
    if not os.path.exists(coeff_data_file):
        #generate fit coefficients
        get_ramp_slope(frame_list, superbias, calFile, mask, slc, degrees, saturation)
    coeff_data = fits.getdata(coeff_data_file)
    x = np.arange(len(frame_list))
    fit_coeff = coeff_data.reshape(2, 4088*4088)
    val_fit = evaluate_poly_array(np.flip(fit_coeff, axis=0), x, poly_type = 'legendre')
    return val_fit.reshape(len(frame_list), 4088, 4088)

def residuals(frame_list, superbias, calFile, mask, slc=((4, 4092), (4, 4092)), degrees=1, saturation=50000):
    residual_file = 'residuals.fits'
    if os.path.exists(residual_file):
        residuals_cube = fits.getdata(residual_file)
    else:
        y_cb = fits.getdata(full_file_list[0].replace('.fits.fz', '.y_cb.fits'))
        y_new = y_cb.reshape(len(frame_list), 4088, 4088)
        residuals_cube = y_new - fit_cube
        fits.writeto(residual_file, residuals_cube, overwrite=True)
    return residuals_cube

fit_cube_file = full_file_list[0].replace('.fits.fz', f'fit_cube_val.fits')
if os.path.exists(fit_cube_file):
    fit_cube = fits.getdata(fit_cube_file)
else:
    fit_cube = val(full_file_list, supercpy, calFile, mask, degrees=3)
    fits.writeto(fit_cube_file, fit_cube, overwrite=True)

res = residuals(full_file_list, supercpy, calFile, mask, degrees=3)
std = np.nanstd(res)
fits.writeto(r'D:\NLC\C1\residuals_cb.fits', res, overwrite=True)

means = []
rms_vals = []
for frame in file_list:
    data = fits.getdata(frame)
    means.append(np.mean(data))
    rms_vals.append(np.sqrt(np.mean(data**2)))

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