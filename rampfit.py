
import numpy as np
from Demos.BackupRead_BackupWrite import outfile
from astropy.io import fits
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
r = (1124972, 1124972+10)

file_list = [file_format.format(n) for n in range(*r)]
full_file_list = file_list
print(len(file_list))

super_bias = fits.getdata('IRRC_calfiles\\super_biasC1.fits.ramp.20231012')
calFile = r'IRRC_calfiles\irrc_weights_C1.h5'

maskFile = fits.getdata(r'IRRC_calfiles\C1_bad_ref_pix_mask.fits')
mask = maskFile>0

supercpy = super_bias.copy()
supercpy[mask[:,:4096]] = 0


def custom_curve_fit(x, y, degrees):
    fitter = fitting.LinearLSQFitter()
    poly = Polynomial1D(degree=degrees)
    return fitter(x,y,poly)

def f(x, m, b):
    return m*x + b

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

# plt.show()
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
def process_files(frame_list,superbias, calFile,multiThread=True, externalPixelFlags:np.ndarray=None):
    for frame in frame_list:
        # frame = full_file_list[i]
        #
        # output_file = frame
        print('number:',frame)

        out_img = irrc_correct_frame_file(frame,superbias,calFile, multiThread, externalPixelFlags)
        outpt = frame + ".out.fits"
        fits.writeto(outpt,out_img,overwrite=True)




def get_ramp_slope(frame_list,superbias, calFile, mask, slc=((4,4092), (4,4092)),degrees=4, saturation=50000):  # ((y1,y2), (x1,x2))
    slopes = []
    corrected_images = []
    y=[]
    for frame in frame_list:
        out_img = irrc_correct_frame_file(frame,superbias,calFile)
        if slc is not None:
            print('slc')
            print(out_img.shape)
            out_img = out_img[:, slc[0][0]:slc[0][1], slc[1][0]:slc[1][1]]
            print(out_img.shape)

        #ISSUE HERE
        saturation_mask = out_img > saturation
        print(saturation_mask.sum()/np.prod(out_img.shape))
        # out_img[saturation_mask] = np.nan

        y.append(out_img)

    x = np.arange(len(y))  #???
    # x = np.asarray([i*np.ones(out_img.shape) for i in x])
    y = np.asarray(y)
    y = y.reshape((x.shape[0],-1))
    sat_pix = (y > saturation).astype(float)
    # x = x.reshape((x.shape[0], -1))
    # output_coeff = np.asarray([np.ones(x.shape) for i in range(degrees+1)]) * np.nan
    # y = y.transpose()
    # x = x.transpose()
    #y = np.expand_dims(saturation_mask,axis=1)

    idx = np.isfinite(y)
    # sum_idx = np.sum(idx,axis=0)
    # for i in np.arange(degrees, sum_idx.max()):
    #     x = np.arange(i)
    #     idy = sum_idx == i
    #     fit_y = np.asarray([frame[idy] for frame in y[:i]])
    #     print('idy', idy.shape)
    #     coefficients, cov_mat = np.polyfit(x, fit_y, degrees, cov=True)
    #     coefficients = np.asarray(coefficients)
    #     print('coefficients', coefficients.shape)
    #     output_coeff[idy] = coefficients

    # fit = np.polyfit(x[idx], y[idx], 1)

    # slope = coefficients[0]   #the slope is just the coefficient of the first degree term (y=mx+b)

    coeff_images = []
    cov_matrices = []
    coefficient_array = []
    cov_array = []
    sum_sat = sat_pix.sum(axis=0)
    bad_pix = sum_sat > 0
    y[:, bad_pix] = 0
    coefficients, cov_mat = np.polyfit(x, y, degrees, cov=True)
    # for i in range(y.shape[1]):
    #     w=unsat_pix[:, i]
    #     if w.sum() < degrees+1:
    #         cov_array.append(np.asarray([np.nan for i in x]))
    #         coefficient_array.append([np.nan for d in range(degrees+1)])
    #     else:
    #         coefficients, cov_mat = np.polyfit(x, y[:, i], degrees,cov=True, w=w)
    #         cov_array.append(cov_mat)
    #         coefficient_array.append(coefficients)
    # coefficient_array = np.asarray(coefficient_array)
    # coefficient_array[bad_pix] = np.asarray([np.nan for d in range(degrees+1)])
    # cov_array = np.asarray(cov_array)
    # cov_array[bad_pix] = np.nan


    # coefficients, cov_mat = curve_fit(f,x[idx],y[idx], nan_policy='omit')
    #coefficients = custom_curve_fit(x[idx], y[idx], degrees)
    coefficients[:, bad_pix] = np.nan

    #ISSUE HERE WITH THE COVARIANCE
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
coeff_images = get_ramp_slope(full_file_list,supercpy,calFile,mask,degrees=1)

#%%
def val(frame_list,superbias, calFile, mask, slc=((4,4092), (4,4092)),degrees=4, saturation=50000):
    for frame in frame_list:
        out_img = irrc_correct_frame_file(frame,superbias,calFile)
        if slc is not None:
            print('slc')
            print(out_img.shape)
            out_img = out_img[:, slc[0][0]:slc[0][1], slc[1][0]:slc[1][1]]
            print(out_img.shape)

        #ISSUE HERE
        saturation_mask = out_img > saturation
        print(saturation_mask.sum()/np.prod(out_img.shape))
        # out_img[saturation_mask] = np.nan

        y.append(out_img)

    x = np.arange(len(y))  #???
    # x = np.asarray([i*np.ones(out_img.shape) for i in x])
    y = np.asarray(y)
    y = y.reshape((x.shape[0],-1))
    for frame in frame_list:
        coeff_data = fits.getdata(full_file_list[0].replace('.fits.fz', f'.img_cb_1deg.fits'))
        x = np.arange(coeff_data.shape[0])
        coefficients = np.polyfit(x, y, degrees, cov=True)
        val = np.polyval(coeff_data,x)
        if slc is not None:
            val = val[:, slc[0][0]:slc[0][1], slc[1][0]:slc[1][1]]
            # val_cube = np.array([fits.getdata(value) for value in val])
            return val

val_cube = val(full_file_list,supercpy,calFile,mask,degrees=1)
residuals_cube = coeff_images - val_cube #this didn't work because shapes are different

fits.writeto('residuals.fits', residuals_cube, overwrite=True)
plt.figure()

summed_residuals = np.sum(residuals_cube, axis=0)

# Save the summed residuals to a new FITS file
fits.writeto('summed_residuals.fits', summed_residuals, overwrite=True)

#ravel flattens pixels into a single array so we can easily put it into a histogram,although im not sure if this is what we want here
plt.hist(summed_residuals.ravel())
plt.title('Histogram of Summed Residuals')
plt.xlabel('Residual Value')
plt.ylabel('Frequency')
plt.grid(True)
plt.show()

#%%

