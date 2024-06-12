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

# open the two images
image1 = fits.getdata(r'D:\NLC\C1\01124973C1_ircc.fits')
image2 = fits.getdata(r'D:\NLC\C1\01124972C1_ircc.fits')

print(image1.shape)
print(image2.shape)

im_new = image1-image2

fits.writeto('subtracted1_image.fits',im_new,overwrite=True)
print(outfile)

file_format = r'D:\NLC\C1\{0:08d}C1.fits.fz'
# r = (1124972, 1125497)
r = (1124972, 1124972+4)

file_list = [file_format.format(n) for n in range(*r)]
full_file_list = file_list
print(len(file_list))

super_bias = fits.getdata('IRRC_calfiles\\super_biasC1.fits.ramp.20231012')
calFile = r'IRRC_calfiles\irrc_weights_C1.h5'

maskFile = fits.getdata(r'IRRC_calfiles\C1_bad_ref_pix_mask.fits')
mask = maskFile>0

supercpy = super_bias.copy()
supercpy[mask[:,:4096]] = 0

#%%

#unused

# def custom_curve_fit(x, y, degrees):
#     fitter = fitting.LinearLSQFitter()
#     poly = Polynomial1D(degree=degrees)
#     return fitter(x,y,poly)

# def f(x, m, b):
#     return m*x + b
#%%

def irrc_correct_frame(
        dataIn, superbias, calFile, multiThread:bool=True, externalPixelFlags:np.ndarray=None,
        superbias_corrected:bool=True
):
    """Returns image cube corrected using IRRC in memory, given calibration file."""

    #Grab coefficients from calFile
    alpha, gamma, zeta = _load_coeffs_from_file(calFile)
    #De-stripe input (subtract off superbias)
    dataIn_tmp = destripe(dataIn[None,:,:])
    dataIn_tmp[:,:,:4096] = dataIn_tmp[:,:,:4096]-superbias[None,:,:]
    #Apply to image cube (after correction add superbias back in)
    dataIRRC = apply_in_memory(dataIn_tmp, alpha, gamma, zeta, multiThread, externalPixelFlags,
                               superbias_corrected)+superbias
    return dataIRRC

def plot_temperatures(frame_list):
    print(fits.getheader(frame_list[0]))
    temps = [fits.open(f)[1].header['TEMPASIC'] for f in frame_list]
    plt.plot(temps)
    # plt.show()

plot_temperatures(file_list)
def irrc_correct_frame_file(filename,superbias, calFile,multiThread:bool=True, externalPixelFlags:np.ndarray=None,
                            superbias_corrected:bool=True):
    dataIn = fits.getdata(filename)
    print(dataIn.shape)
    dataIn = dataIn[:, 6:]
    dcpy = dataIn.copy()
    dcpy[mask]=0
    print(dataIn.shape)
    return irrc_correct_frame(dcpy, superbias, calFile,multiThread, externalPixelFlags, superbias_corrected)

print(irrc_correct_frame_file(r'D:\NLC\C1\01124973C1.fits.fz',supercpy,calFile))
def get_ramp_cds(frame_list,superbias, calFile,multiThread=True, externalPixelFlags:np.ndarray=None,
                 superbias_corrected:bool=True):
    image1 = irrc_correct_frame_file(frame_list[0], superbias, calFile, multiThread, externalPixelFlags, superbias_corrected)
    image2 = irrc_correct_frame_file(frame_list[-1], superbias, calFile, multiThread, externalPixelFlags, superbias_corrected)
    cds = image2-image1
    f = frame_list[0]+".cds.fits"
    fits.writeto(f, cds, overwrite=True)
    return f
#%%
def get_ramp_slope(frame_list,superbias, calFile, mask, slc=((4,4092), (4,4092)),degrees=1, saturation=50000):  # ((y1,y2), (x1,x2))
    y=[]
    for frame in frame_list:
        out_img = irrc_correct_frame_file(frame,superbias,calFile)
        if slc is not None:
            print('slc')
            print(out_img.shape)
            out_img = out_img[:, slc[0][0]:slc[0][1], slc[1][0]:slc[1][1]]
            print(out_img.shape)

        saturation_mask = out_img > saturation
        print(saturation_mask.sum()/np.prod(out_img.shape))

        y.append(out_img)

    x = np.arange(len(y))  #???
    # x = np.asarray([i*np.ones(out_img.shape) for i in x])
    y = np.asarray(y)
    fits.writeto(frame_list[0].replace('.fits.fz', f'.y_cb.fits'), y, overwrite=True)
    y = y.reshape((x.shape[0],-1))
    sat_pix = (y > saturation).astype(float)

    coeff_images = []
    cov_matrices = []
    sum_sat = sat_pix.sum(axis=0)
    bad_pix = sum_sat > 0
    y[:, bad_pix] = 0
    #coefficients, cov_mat = np.polyfit(x, y, degrees, cov=True)

    coefficients,cov_mat = np.polynomial.legendre.legfit(x,y,degrees,cov=True)

    coefficients[:, bad_pix] = np.nan

    for coeff in coefficients:
        coeff_image = coeff.reshape(out_img.shape)
        coeff_images.append(coeff_image)
        cov_mat_img = coeff.reshape(out_img.shape)
        cov_matrices.append(cov_mat_img)

    coeff_images = np.stack(coeff_images,axis=0)
    cov_matrices = np.stack(cov_matrices,axis=0)

    fits.HDUList([fits.PrimaryHDU(coeff_images)]).writeto(frame_list[0].replace('.fits.fz', f'.img_cb_{degrees}deg.fits'), overwrite=True)
    fits.HDUList([fits.PrimaryHDU(cov_matrices)]).writeto(frame_list[0].replace('.fits.fz', f'.pcov_cb_{degrees}deg.fits'), overwrite=True)

    return coeff_images, cov_matrices

print('get_ramp_slope')
#coeff_images = get_ramp_slope(full_file_list,supercpy,calFile,mask,degrees=1)

#%%
def evaluate_poly_array(coeffs,x_array):
    output_arrays = []
    for x in x_array:
        output_array = np.zeros(coeffs.shape[1])
        for n, coeff in enumerate(coeffs):
            x_array = x**n
            output_array += coeff*x_array
        output_arrays.append(output_array)
    return np.asarray(output_arrays)

def val(frame_list, superbias, calFile, mask, slc=((4, 4092), (4, 4092)), degrees=1, saturation=50000):
    coeff_data = fits.getdata(full_file_list[0].replace('.fits.fz', f'.img_cb_{degrees}deg.fits'))

    x = np.arange(len(frame_list))

    fit_coeff = coeff_data.reshape(2,4088*4088)

    val_fit = evaluate_poly_array(np.flip(fit_coeff,axis=0),x)

    return val_fit.reshape(len(frame_list), 4088, 4088)

fit_cube = val(full_file_list, supercpy, calFile, mask, degrees=1)

values = fits.writeto(full_file_list[0].replace('.fits.fz', f'fit_cube_val.fits'), fit_cube, overwrite=True)

fit_cube = val(full_file_list, supercpy, calFile, mask, degrees=1)
#%%
def residuals(frame_list, superbias, calFile, mask, slc=((4, 4092), (4, 4092)), degrees=1, saturation=50000):
    y_cb = fits.getdata(full_file_list[0].replace('.fits.fz', '.y_cb.fits'))

    y_new = y_cb.reshape(len(frame_list), 4088, 4088)

    residuals_cube = y_new[0] - fit_cube[0]

    fits.writeto('residuals.fits', residuals_cube, overwrite=True)

    return residuals_cube
#%%
print('PRINTED')
res = residuals(full_file_list, supercpy, calFile, mask, degrees=1)
print(np.nanstd(res))
# summed_residuals = np.sum(res, axis=0)
#%%
res_cb = fits.writeto('residuals_cb.fits', res,overwrite=True)

#summed resid for all, aka last frame
plt.figure()
hist = np.histogram(res[np.isfinite(res)],bins=2)

plt.plot(hist[0],hist[1],color='blue')
plt.title('Histogram of Residuals')
plt.xlabel('Residual Value')
plt.ylabel('Frequency')
plt.grid(True)
plt.show()

# #last frame
# last_frame = res[-1]
# plt.figure()
# plt.hist(last_frame.ravel(),bins=50,color='red')
# plt.xlim(-200,200)
# plt.title('Histogram of Last Frame')
# plt.xlabel('Residual Value')
# plt.ylabel('Frequency')
# plt.grid(True)
# plt.show()

#%%
