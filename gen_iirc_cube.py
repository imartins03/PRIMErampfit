from astropy.io import fits
import numpy as np
from irrc.apply_correction import _load_coeffs_from_file, apply_in_memory
from irrc.util import destripe
import matplotlib.pyplot as plt

file_format = r'D:\NLC\C1\{0:08d}C1.fits.fz'
r = (1124972, 1124972+100)
file_list = [file_format.format(n) for n in range(*r)]  # Generate list of file paths using the specified range
full_file_list = file_list
super_bias_path = 'IRRC_calfiles/super_biasC1.fits.ramp.20231012'
calFile = r'IRRC_calfiles\irrc_weights_C1.h5'
maskFile_path = r'IRRC_calfiles/C1_bad_ref_pix_mask.fits'
y_cube_path = r'D:\NLC\C1\y_cube_100.fits'

super_bias = fits.getdata(super_bias_path)  # Load super bias data
maskFile = fits.getdata(maskFile_path)  # Load mask file data
mask = maskFile > 0  # Create mask for bad pixels
supercpy = super_bias.copy()  # Create a copy of the super bias
supercpy[mask[:, :4096]] = 0  # Apply mask to the super bias copy

def irrc_correct_frame(dataIn, superbias, calFile, multiThread=True, externalPixelFlags=None, superbias_corrected=True):
    """Returns image cube corrected using IRRC in memory, given calibration file."""
    alpha, gamma, zeta = _load_coeffs_from_file(calFile)  # Load IRRC coefficients
    dataIn_tmp = destripe(dataIn[None, :, :])  # Destripe the input data
    dataIn_tmp[:, :, :4096] -= superbias[None, :, :]  # Subtract the super bias
    dataIRRC = apply_in_memory(dataIn_tmp, alpha, gamma, zeta, multiThread, externalPixelFlags, superbias_corrected) + superbias  # Apply IRRC correction
    return dataIRRC

def irrc_correct_frame_file(filename, superbias, calFile, multiThread=True, externalPixelFlags=None, superbias_corrected=True):
    dataIn = fits.getdata(filename)  # Load data from FITS file
    dataIn = dataIn[:, 6:]  # Adjust data by removing first 6 columns?
    dcpy = dataIn.copy()  # Create a copy of the data
    dcpy[mask] = 0  # Apply mask to the data copy
    return irrc_correct_frame(dcpy, superbias, calFile, multiThread, externalPixelFlags, superbias_corrected)  # Correct the frame using IRRC

def generate_y_cube(frame_list, superbias, calFile, slc=((4, 4092), (4, 4092))):
    y = []  # Initialize list to hold corrected frames

    for frame in frame_list:  # Loop through each frame in the list
        out_img = irrc_correct_frame_file(frame, superbias, calFile)  # Correct frame using IRRC and calibration file
        if slc is not None:
            out_img = out_img[:, slc[0][0]:slc[0][1], slc[1][0]:slc[1][1]]  # Slice the corrected frame
        y.append(out_img)  # Append corrected frame to list
    y = np.asarray(y)  # Convert list to NumPy array
    fits.writeto(y_cube_path, y, overwrite=True)  # Write the array to a FITS file
    return y

y_cube = generate_y_cube(full_file_list, supercpy, calFile)  # Generate the y cube
