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

fits.writeto('IRRC_calfiles/subtracted1_image.fits', im_new, overwrite=True)
print(outfile)
#%%

file_format = r'D:\NLC\C1\{0:08d}C1.fits.fz'
# r = (1124972, 1125497)
r = (1124972, 1124972+10)
file_list = [file_format.format(n) for n in range(*r)]
full_file_list = file_list
print(file_list)
super_bias = fits.getdata(r'C:\\PycharmProjects\\PRIMErampfit\\IRRC_calfiles\\super_biasC1.fits.ramp.20231012')
calFile = r'C:\PycharmProjects\PRIMErampfit\IRRC_calfiles\irrc_weights_C1.h5'

maskFile = fits.getdata(r'C:\\PycharmProjects\\PRIMErampfit\\IRRC_calfiles\\C1_bad_ref_pix_mask.fits')
bad_pixel_mask = maskFile>0

supercpy = super_bias.copy()
supercpy[bad_pixel_mask[:,:4096]] = 0

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
    dcpy[externalPixelFlags]=0
    print(dataIn.shape)
    return irrc_correct_frame(dcpy, superbias, calFile,multiThread, externalPixelFlags, superbias_corrected)

plt.show()
print(irrc_correct_frame_file(r'D:\NLC\C1\01124973C1.fits.fz',supercpy,calFile))
#%%
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
def get_ramp_slope(frame_list,superbias, calFile, externalPixelFlags, slc=((4,4088), (4,4088)),degrees=8):  # ((y1,y2), (x1,x2))
    slopes = []
    corrected_images = []
    y=[]
    # filename, superbias, calFile, multiThread: bool = True, externalPixelFlags:np.ndarray = None,
    # superbias_corrected: bool = True):
    for frame in frame_list:
        out_img = irrc_correct_frame_file(frame,superbias,calFile,externalPixelFlags=externalPixelFlags, superbias_corrected=True)
        if slc is not None:
            print('slc')
            print(out_img.shape)
            out_img = out_img[:, slc[0][0]:slc[0][1], slc[1][0]:slc[1][1]]
            print(out_img.shape)
        y.append(out_img)
    x = np.arange(len(y))  #???
    y = np.asarray(y)

    fits.writeto(frame_list[0].replace('.fits.fz', f'.y_cb.fits'), y, overwrite=True)

    y = y.reshape((x.shape[0],-1))
    print(y.ndim)
    # print('the output image array is:',out_img)

    coeff_images = []


    coefficients = np.polyfit(x,y,degrees)
    coefficients = np.flip(coefficients,axis=0)
    for coeff in coefficients:
        coeff_image = coeff.reshape(out_img.shape)
        coeff_images.append(coeff_image)

    coeff_images = np.stack(coeff_images,axis=0)

    fits.writeto(frame_list[0].replace('.fits.fz', f'.img_cb_{degrees}deg.fits'), coeff_images, overwrite=True)  ##
    # fits.writeto(frame_list[0].replace('.fits.fz', '.ramp.fits'), fit_image,overwrite=True)
    return coeff_images

print('get_ramp_slope')
coeff_images = get_ramp_slope(full_file_list,supercpy,calFile,bad_pixel_mask,degrees=1)
# print(coeff_images)


# def get_ramp_slope(frame_list,superbias, calFile, mask, slc=((4,4088), (4,4088)),degrees=8):
#     slopes = []
#     corrected_images = []
#     y=[]
#     for frame in frame_list:
#         out_img = irrc_correct_frame_file(frame,superbias,calFile)
#         y.append(out_img)
#
#     x = np.arange(len(y))  # ???
#     y = np.asarray(y)
#
#     fits.writeto(frame_list[0].replace('.fits.fz', f'.y_cb.fits'), y, overwrite=True)
#
#     y = y.reshape((x.shape[0], -1))
#     print(y.ndim)
#     # print('the output image array is:',out_img)
#
#     coeff_images = []
#
#     coefficients = np.polyfit(x, y, degrees)
#     for coeff in coefficients:
#         coeff_image = coeff.reshape(out_img.shape)
#         coeff_images.append(coeff_image)
#
#     coeff_images = np.stack(coeff_images, axis=0)
#
#     fits.writeto(frame_list[0].replace('.fits.fz', f'.img_cb_{degrees}deg.fits'), coeff_images, overwrite=True)  ##
#     # fits.writeto(frame_list[0].replace('.fits.fz', '.ramp.fits'), fit_image,overwrite=True)
#     return coeff_images
#
# print('get_ramp_slope')
# coeff_images = get_ramp_slope(full_file_list, supercpy, calFile, mask, degrees=1)
#%%

import numpy as np
import matplotlib.pyplot as plt

# Assuming coeff_images is a numpy array where each coefficient is stored as a separate image
# Assuming coeff_images.shape = (n_coefficients, height, width)

# Extracting the first order coefficients
first_order_coefficients = coeff_images[0]  # Assuming the first coefficient is the first order

# Number of orders (excluding the first order)
n_orders = coeff_images.shape[0] - 1

# Create a 2x3 grid of subplots
fig, axs = plt.subplots(2, 3, figsize=(12, 8), sharex=True)

# Plot each order coefficient vs. the first order in a separate subplot
for order, ax in zip(range(1, min(n_orders + 1, 8)), axs.flatten()):
    order_coefficients = coeff_images[order]  # Assuming coefficients for each order are stored sequentially
    # Reshape the coefficients to a 1D array (flatten the image)
    order_coefficients_flat = order_coefficients.flatten()
    # Plot the coefficients vs. the first order
    ax.scatter(first_order_coefficients.flatten(), order_coefficients_flat, s=1)  # Adjust marker size as needed
    ax.set_ylabel(f'Order {order} Coefficient')
    ax.set_title(f'Order {order}')

# Add labels and title to the x-axis (assuming all subplots share the same x-axis)
fig.text(0.5, 0.04, 'First Order Coefficient', ha='center', va='center')
fig.text(0.04, 0.5, 'Order Coefficient', ha='center', va='center', rotation='vertical')

# Adjust layout to prevent overlap
plt.tight_layout()

# Show the plot
plt.show()

