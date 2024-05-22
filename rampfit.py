
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



# open the two images
image1 = fits.getdata(r'D:\NLC\C1\01124973C1_ircc.fits')
image2 = fits.getdata(r'D:\NLC\C1\01124972C1_ircc.fits')

print(image1.shape)
print(image2.shape)

im_new = image1-image2

fits.writeto('subtracted1_image.fits',im_new,overwrite=True)
print(outfile)
#%%

image1 = fits.getdata(r'D:\NLC\C1\01124973C1_ircc.fits')
image2 = fits.getdata(r'D:\NLC\C1\01124972C1_ircc.fits')
image3 = fits.getdata(r'D:\NLC\C1\01124974C1_ircc.fits')
image4 = fits.getdata(r'D:\NLC\C1\01124975C1_ircc.fits')

frame_list = [image1,image2,image3,image4]
superbias = fits.getdata('C:\\PycharmProjects\\PRIMErampfit\\IRRC_calfiles\\super_biasC1.fits.ramp.20231012')
calFile = r'C:\PycharmProjects\PRIMErampfit\IRRC_calfiles\irrc_weights_C1.h5'

img_cube_list = []
hdr_list = []
for f in frame_list:
    hdu = fits.open(f)
    d = hdu[image_hdu].data[:, 6:]  # Assumes data is in 1st HDU of FITS and headers are present
    h = hdu[image_hdu].header
    hdr_list.append(h)
    # Mask bad reference pixels
    dcpy = d.copy()
    dcpy[mask] = 0
    img_cube_list.append(dcpy)
dataIn = np.asarray(img_cube_list)

def irrc_correct_frame(image, superbias, calFile,multiThread:bool=True, externalPixelFlags:np.ndarray=None, superbias_corrected:bool=True):
    """Returns image cube corrected using IRRC in memory, given calibration file."""

    #Grab coefficients from calFile
    alpha, gamma, zeta = _load_coeffs_from_file(calFile)
    #De-stripe input (subtract off superbias)
    dataIn_tmp = destripe(dataIn[None,:,:])
    dataIn_tmp[:,:,:4096] = dataIn_tmp[:,:,:4096]-superbias[None,:,:]
    #Apply to image cube (after correction add superbias back in)
    dataIRRC = apply_in_memory(dataIn_tmp, alpha, gamma, zeta, multiThread, externalPixelFlags, superbias_corrected=True)+superbias
    return dataIRRC

def get_ramp_cds(frame_list,superbias, calFile,multiThread=True, externalPixelFlags:np.ndarray=None, superbias_corrected:bool=True):
    cds = irrc_correct_frame((frame_list[-1]-frame_list[0]),superbias,calFile, multiThread, externalPixelFlags, superbias_corrected)
    fits.writeto('cds_file.fits', cds, overwrite=True)

corrected_frame = irrc_correct_frame(dataIn, superbias, calFile)
fits.writeto('corrected_frame_image1.fits', corrected_frame[0], overwrite=True)



#%%
print(irrc_correct_frame(image1,superbias,calFile))



#def (image1,)
 #   np.linearregression

#def (image1,)
 #   return


#%%

#Number of bits in the ADC
#p = 8

#z is the max number of samples
#z = 1894
#q = math.log2(z)





