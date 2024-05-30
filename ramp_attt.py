import numpy as np
from astropy.io import fits
from irrc.apply_correction import _load_coeffs_from_file, apply_in_memory
from irrc.util import destripe

file_format = r'D:\NLC\C1\{0:08d}C1.fits.fz'
# r = (1124972, 1125497)
r = (1124972, 1124972+2)
file_list = [file_format.format(n) for n in range(*r)]
full_file_list = file_list
print(file_list)
super_bias = fits.getdata('C:\\PycharmProjects\\PRIMErampfit\\IRRC_calfiles\\super_biasC1.fits.ramp.20231012')
calFile = r'C:\PycharmProjects\PRIMErampfit\IRRC_calfiles\irrc_weights_C1.h5'

maskFile = fits.getdata(r'C:\PycharmProjects\PRIMErampfit\IRRC_calfiles\C1_bad_ref_pix_mask.fits')
mask = maskFile>0

supercpy = super_bias.copy()
supercpy[mask[:,:4096]] = 0

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
    maskFile = fits.getdata('C:\\PycharmProjects\\PRIMErampfit\\IRRC_calfiles\\C1_bad_ref_pix_mask.fits')
    mask = maskFile > 0
    dcpy = dataIn.copy()
    dcpy[mask] = 0
    return irrc_correct_frame(dcpy, superbias, calFile, multiThread, externalPixelFlags, superbias_corrected)

def get_ramp_slope(frame_list, superbias, calFile, mask, degrees=8, slc=((4, 4088), (4, 4088))):
    y = []
    for frame in frame_list:
        out_img = irrc_correct_frame_file(frame, superbias, calFile)
        if slc is not None:
            out_img = out_img[:, slc[0][0]:slc[0][1], slc[1][0]:slc[1][1]]
        y.append(out_img)
    x = np.arange(len(y))
    y = np.asarray(y)
    y = y.reshape((x.shape[0], -1))

    coeffs = []
    for degree in range(1, degrees + 1):
        coefficients = np.polyfit(x, y, degree)
        coeffs.append(coefficients)

    coeffs = np.asarray(coeffs)
    coeff_images = []
    for i in range(coeffs.shape[0]):
        for j in range(coeffs.shape[1]):
            coeff_image = coeffs[i, j].reshape(y.shape[1:])
            fits.writeto(frame_list[0].replace('.fits.fz', f'.coeff_{i}_{j}.fits'), coeff_image, overwrite=True)
            coeff_images.append(coeff_image)

    return coeff_images

degrees = 8
coeff_images = get_ramp_slope(full_file_list,supercpy,calFile,mask,degrees=2)
