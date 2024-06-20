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


y_cube_path = r'D:\NLC\C1\y_cube.fits'
fit_cube_path = r'D:\NLC\C1\fit_cube.fits'
residuals_cube_path = r'D:\NLC\C1\residuals.fits'

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

y_cube = generate_y_cube(full_file_list, supercpy, calFile)  #call it once with this in the y function


#