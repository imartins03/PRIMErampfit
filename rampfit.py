
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



# open the two images
image1 = fits.getdata(r'D:\NLC\C1\01124973C1_ircc.fits')
image2 = fits.getdata(r'D:\NLC\C1\01124972C1_ircc.fits')

print(image1.shape)
print(image2.shape)

im_new = image1-image2

fits.writeto('subtracted1_image.fits',im_new,overwrite=True)
print(outfile)
#%%

# image1 = fits.getdata(r'D:\NLC\C1\01124973C1_ircc.fits')
# image2 = fits.getdata(r'D:\NLC\C1\01124972C1_ircc.fits')
# image3 = fits.getdata(r'D:\NLC\C1\01124974C1_ircc.fits')
# image4 = fits.getdata(r'D:\NLC\C1\01124975C1_ircc.fits')
#
# frame_list = [image1,image2,image3,image4]
file_format = r'D:\NLC\C1\{0:08d}C1.fits.fz'
# r = (1124972, 1125497)
r = (1124972, 1124972+10)
file_list = [file_format.format(n) for n in range(*r)]
full_file_list = file_list
print(file_list)
super_bias = fits.getdata('C:\\PycharmProjects\\PRIMErampfit\\IRRC_calfiles\\super_biasC1.fits.ramp.20231012')
calFile = r'C:\PycharmProjects\PRIMErampfit\IRRC_calfiles\irrc_weights_C1.h5'

maskFile = fits.getdata(r'C:\PycharmProjects\PRIMErampfit\IRRC_calfiles\C1_bad_ref_pix_mask.fits')
mask = maskFile>0

supercpy = super_bias.copy()
supercpy[mask[:,:4096]] = 0



# img_cube_list = []
# hdr_list = []
# for f in frame_list:
#     hdu = fits.open(f)
#     d = hdu[image_hdu].data[:, 6:]  # Assumes data is in 1st HDU of FITS and headers are present
#     h = hdu[image_hdu].header
#     hdr_list.append(h)
#     # Mask bad reference pixels
#     dcpy = d.copy()
#     dcpy[mask] = 0
#     img_cube_list.append(dcpy)
# dataIn = np.asarray(img_cube_list)


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
#%%
# file_list = [r'D:\NLC\C1\01124972C1.fits.fz',r'D:\NLC\C1\01124973C1.fits.fz']
#,r'D:\NLC\C1\01124974C1.fits.fz',r'D:\NLC\C1\01124975C1.fits.fz']

def get_ramp_cds(frame_list,superbias, calFile,multiThread=True, externalPixelFlags:np.ndarray=None,
                 superbias_corrected:bool=True):
    image1 = irrc_correct_frame_file(frame_list[0], superbias, calFile, multiThread, externalPixelFlags, superbias_corrected)
    image2 = irrc_correct_frame_file(frame_list[-1], superbias, calFile, multiThread, externalPixelFlags, superbias_corrected)
    cds = image2-image1
    f = frame_list[0]+".cds.fits"
    fits.writeto(f, cds, overwrite=True)
    return f

# print(get_ramp_cds(file_list,supercpy,calFile,superbias_corrected=True))
#get_ramp_cds(file_list,supercpy,calFile)


# corrected_frame = irrc_correct_frame(dataIn, superbias, calFile)
# fits.writeto('corrected_frame_image1.fits', corrected_frame[0], overwrite=True)

#%%

#Number of bits in the ADC
#p = 8

#z is the max number of samples
#z = 1894
#q = math.log2(z)

#make a new function but ramp reduced
#same inputs and write a file at the end, but instead of subtraction between teo frames do a linear fit of the frames, use linear regression numpy function

#make a loop to call out/generate images, create a linear fit, and then use np.poly fit (with the 1 so its a 1 degreee polynomial) to create fit)

# full_file_list = [r'D:\NLC\C1\01124972C1.fits.fz',r'D:\NLC\C1\01124973C1.fits.fz',r'D:\NLC\C1\01124974C1.fits.fz',r'D:\NLC\C1\01124975C1.fits.fz']

# def process_files(frame_list,superbias, calFile,multiThread=True, externalPixelFlags:np.ndarray=None):
#     for i in range(len(full_file_list)):
#         frame = full_file_list[i]
#
#         output_file = frame
#
#         cor_img = irrc_correct_frame_file(output_file,superbias,calFile,mask)
#         outpt = frame + ".cor.fits"
#         fits.writeto(outpt,cor_img,overwrite=True)
#         print(f'The output is saved to {outpt}')

def process_files(frame_list,superbias, calFile,multiThread=True, externalPixelFlags:np.ndarray=None):
    for frame in frame_list:
        # frame = full_file_list[i]
        #
        # output_file = frame
        print('number:',frame)

        out_img = irrc_correct_frame_file(frame,superbias,calFile, multiThread, externalPixelFlags)
        outpt = frame + ".out.fits"
        fits.writeto(outpt,out_img,overwrite=True)
        # print(f'The output is saved to {outpt}')

# print(frame_list)
# print(process_files(full_file_list,supercpy,calFile))

#%%

def get_ramp_slope(frame_list,superbias, calFile, mask, slc=((4,4088), (4,4088))):  # ((y1,y2), (x1,x2))
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
        y.append(out_img)
    x = np.arange(len(y))  #???
    y = np.asarray(y)
    y = y.reshape((x.shape[0],-1))
    coefficients = np.polyfit(x,y,2)
    slope = coefficients[0]   #the slope is just the coefficient of the first degree term (y=mx+b)
    # slopes.append(slope)   #adding the slopes to the list

    # #To create a linear fit line, polyval evaluates specific values for coefficients, polyfit fits and return the actual coefficients
    # fit = np.polyval(coefficients,x)
    fit_image = slope.reshape(out_img.shape)  #tried to fit it back to the image here (instead of flattened)
    print(fit_image)
    # corrected_images.append(cor_img)
    fits.writeto(frame_list[0].replace('.fits.fz', '.ramp.fits'), fit_image,overwrite=True)
    return slopes,corrected_images
print('get_ramp_slope')
slopes,corrected_images = get_ramp_slope(full_file_list,supercpy,calFile,mask)
# print(f'slopes={fit_image},cor_imgs={corrected_images}')

#%%

# #NEED TO MAKE SURE ALL IMAGES DOWNLOAD AND NOT JUST ONE
# def corrected_files(frame_list,superbias, calFile,multiThread=True, externalPixelFlags:np.ndarray=None):
#     slopes, corrected_images = get_ramp_slope(frame_list, superbias, calFile, mask)
#     for frame,cor_img in zip(frame_list,corrected_images):
#         print('number processing:',frame)
#
#         # img = get_ramp_slope(frame_list,superbias, calFile, mask)
#         output = frame + ".cor.fits"
#         fits.writeto(output, cor_img, overwrite=True)
#
# print(corrected_files(full_file_list,supercpy,calFile))

#%%



