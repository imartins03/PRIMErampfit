
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

super_bias = fits.getdata('C:\\PycharmProjects\\PRIMErampfit\\IRRC_calfiles\\super_biasC1.fits.ramp.20231012')
calFile = r'C:\PycharmProjects\PRIMErampfit\IRRC_calfiles\irrc_weights_C1.h5'

maskFile = fits.getdata(r'C:\PycharmProjects\PRIMErampfit\IRRC_calfiles\C1_bad_ref_pix_mask.fits')
mask = maskFile>0

supercpy = super_bias.copy()
supercpy[mask[:,:4096]] = 0

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

def irrc_correct_frame_file(filename,superbias, calFile,multiThread:bool=True, externalPixelFlags:np.ndarray=None,
                            superbias_corrected:bool=True):
    dataIn = fits.getdata(filename)
    print(dataIn.shape)
    dataIn = dataIn[:, 6:]
    dcpy = dataIn.copy()
    dcpy[mask]=0
    print(dataIn.shape)
    return irrc_correct_frame(dcpy, superbias, calFile,multiThread, externalPixelFlags, superbias_corrected)

plt.show()
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

def get_ramp_slope(frame_list,superbias, calFile, mask, slc=((4,4088), (4,4088)),degrees=4):  # ((y1,y2), (x1,x2))
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
    saturation_mask = out_img > 20000
    out_img[saturation_mask] = np.nan

    y.append(out_img)
    x = np.arange(len(y))  #???
    y = np.asarray(y)
    y = y.reshape((x.shape[0],-1))

    # slope = coefficients[0]   #the slope is just the coefficient of the first degree term (y=mx+b)

    coeff_images = []
    cov_matrices = []
    coefficients, cov_mat = np.polyfit(x, y, degrees,cov=True)

    #ISSUE HERE WITH THE COVARIANCE
    for coeff in coefficients:
        coeff_image = coeff.reshape(out_img.shape)
        coeff_images.append(coeff_image)
        cov_mat_img = coeff.reshape(out_img.shape)
        cov_matrices.append(cov_mat_img)

    coeff_images = np.stack(coeff_images,axis=0)
    cov_matrices = np.stack(cov_matrices,axis=0)
    # hdu = fits.ImageHDU(coeff_images)
    fits.HDUList([fits.PrimaryHDU(coeff_images)]).writeto(frame_list[0].replace('.fits.fz', f'.img_cb_{degrees}deg.fits'), overwrite=True)
    fits.HDUList([fits.PrimaryHDU(cov_matrices)]).writeto(frame_list[0].replace('.fits.fz', f'.pcov_cb_{degrees}deg.fits'), overwrite=True)

    # fits.writeto(frame_list[0].replace('.fits.fz', f'.img_cb_{degrees}deg.fits'), coeff_images, overwrite=True)  ##
    # fits.writeto(frame_list[0].replace('.fits.fz', '.ramp.fits'), fit_image,overwrite=True)

    return coeff_images, cov_matrices

print('get_ramp_slope')
coeff_images = get_ramp_slope(full_file_list,supercpy,calFile,mask,degrees=1)



