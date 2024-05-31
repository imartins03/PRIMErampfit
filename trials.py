import numpy as np
from astropy.io import fits
from irrc.apply_correction import _load_coeffs_from_file, apply_in_memory
from irrc.util import destripe
import matplotlib.pyplot as plt

# Load images
image1 = fits.getdata(r'D:\NLC\C1\01124973C1_ircc.fits')

image2 = fits.getdata(r'D:\NLC\C1\01124972C1_ircc.fits')

im_new = image1 - image2
fits.writeto('subtracted1_image.fits', im_new, overwrite=True)

file_format = r'D:\NLC\C1\{0:08d}C1.fits.fz'
r = (1124972, 1124972 + 4)
file_list = [file_format.format(n) for n in range(*r)]
full_file_list = file_list

super_bias = fits.getdata('C:\\PycharmProjects\\PRIMErampfit\\IRRC_calfiles\\super_biasC1.fits.ramp.20231012')

calFile = r'C:\PycharmProjects\PRIMErampfit\IRRC_calfiles\irrc_weights_C1.h5'
maskFile = fits.getdata(r'C:\PycharmProjects\PRIMErampfit\IRRC_calfiles\C1_bad_ref_pix_mask.fits')

mask = maskFile > 0

supercpy = super_bias.copy()
supercpy[mask[:, :4096]] = 0


def irrc_correct_frame(dataIn, superbias, calFile, multiThread=True, externalPixelFlags=None, superbias_corrected=True):
    alpha, gamma, zeta = _load_coeffs_from_file(
        calFile)

    dataIn_tmp = destripe(dataIn[None, :, :])
    dataIn_tmp[:, :, :4096] -= superbias[None, :, :]
    dataIRRC = apply_in_memory(dataIn_tmp, alpha, gamma, zeta, multiThread, externalPixelFlags,
                               superbias_corrected) + superbias
    return dataIRRC


def irrc_correct_frame_file(
        filename, superbias, calFile, multiThread=True, externalPixelFlags=None, superbias_corrected=True):
    dataIn = fits.getdata(filename)
    dataIn = dataIn[:, 6:]
    dcpy = dataIn.copy()
    dcpy[mask] = 0
    return irrc_correct_frame(dcpy, superbias, calFile, multiThread, externalPixelFlags, superbias_corrected)


def get_ramp_slope(frame_list, superbias, calFile, mask, slc=((4, 4088), (4, 4088)), degrees=4):
    corrected_images = []
    y = []

    for frame in frame_list:
        out_img = irrc_correct_frame_file(frame, superbias, calFile)
        if slc is not None:
            out_img = out_img[:, slc[0][0]:slc[0][1], slc[1][0]:slc[1][1]]

        # Mask out pixels above 20,000 and replace them with np.nan
        out_img[out_img > 20000] = np.nan
        y.append(out_img)

    x = np.arange(len(frame_list))  # Adjusted to match the length of frame_list

    y = np.asarray(y)
    y = y.reshape((x.shape[0], -1))

    # Masking NaNs for polynomial fitting
    mask = np.isnan(y)
    y_masked = np.ma.masked_array(y, mask)

    coeff_images = []
    try:
        coefficients, cov_mat = np.polyfit(x, y_masked, degrees, cov=True)
    except np.linalg.LinAlgError as e:
        print("An error occurred during polynomial fitting: ", e)
        return None, None

    for coeff in coefficients:
        coeff_image = coeff.reshape(out_img.shape)
        coeff_images.append(coeff_image)

    coeff_images = np.stack(coeff_images, axis=0)
    fits.HDUList([fits.PrimaryHDU(coeff_images)]).writeto(
        frame_list[0].replace('.fits.fz', f'.img_cb_{degrees}deg.fits'), overwrite=True)

    cov_matrices = []
    for i in range(degrees + 1):
        for j in range(degrees + 1):
            cov_mat_img = cov_mat[i, j].reshape(out_img.shape)
            cov_matrices.append(cov_mat_img)

    cov_matrices = np.stack(cov_matrices, axis=0)
    fits.HDUList([fits.PrimaryHDU(cov_matrices)]).writeto(
        frame_list[0].replace('.fits.fz', f'.cov_cb_{degrees}deg.fits'), overwrite=True)

    return coeff_images, cov_matrices


print('Generating ramp slopes...')
coeff_images, cov_matrices = get_ramp_slope(full_file_list, supercpy, calFile, mask, degrees=4)