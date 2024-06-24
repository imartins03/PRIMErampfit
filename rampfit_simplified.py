from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

super_bias_path = 'IRRC_calfiles\\super_biasC1.fits.ramp.20231012'
calFile = r'IRRC_calfiles\irrc_weights_C1.h5'
maskFile_path = r'IRRC_calfiles\C1_bad_ref_pix_mask.fits'
y_cube_path = r'D:\NLC\C1\y_cube.fits'
fit_cube_path = r'D:\NLC\C1\fit_cube.fits'
fit_coeff_path = r'D:\NLC\C1\fit_coeff.fits'
residuals_cube_path = r'D:\NLC\C1\residuals.fits'

super_bias = fits.getdata(super_bias_path)
maskFile = fits.getdata(maskFile_path)
mask = maskFile > 0
supercpy = super_bias.copy()
supercpy[mask[:, :4096]] = 0

def evaluate_poly_array(coeffs, a_array, poly_type='power'):
    output_arrays = []
    for a in a_array:
        print(a)
        if poly_type == 'power':
            output_array = np.zeros(coeffs.shape[1])
            for n, coeff in enumerate(coeffs):
                output_array += coeff*(a**n)
            output_arrays.append(output_array)
    return np.asarray(output_arrays)

def generate_fit_cube(frame_num, degrees, saturation=50000):
    y_cube = fits.getdata(y_cube_path)
    print(np.ndim(y_cube))
    x = y_cube.shape[0]   #x is the dimension of the data cube (number of frames)
    # x = np.arange(len(file_list))

    y = y_cube.reshape(x, 4088, 4088)
    print(y.shape)
    y = y_cube.reshape(x, -1)

    print(y.shape)

    z = np.arange(len(y))
    print(z)
    # coefficients, _ = np.polynomial.legendre.legfit(z, y, degrees)
    coefficients,_ = np.polyfit(z, y, degrees,cov=True)

    print(coefficients.shape)
    fit_coeff = coefficients.reshape(degrees+1,4088,4088)
    fits.writeto(fit_coeff_path, fit_coeff, overwrite=True)
    print(fit_coeff.shape)
    fit_cube = evaluate_poly_array(np.flip(coefficients, axis=0), z)

    print(y.shape[1])
    print(fit_cube.shape)

    fit_cube = fit_cube.reshape(len(z), 4088, 4088)
    fits.writeto(fit_cube_path, fit_cube, overwrite=True)
    return fit_cube

saturation=50000
degrees = 1
generate_fit_cube(np.linspace(1,4,4), degrees, saturation)  #HAVE TO CHANGE

def calculate_residuals():
    y_cube = fits.getdata(y_cube_path)
    y_cube = y_cube[:,0,:,:]
    fit_cube = fits.getdata(fit_cube_path)
    residuals_cube = y_cube - fit_cube
    fits.writeto(residuals_cube_path, residuals_cube, overwrite=True)
    return residuals_cube

res = calculate_residuals()

# frame_num = []
# means = []
# rms_vals = []
# residual_vals = []
#
# for i in range(res.shape[0]):
#     data = res[i]
#     means.append(np.mean(data))
#     rms_vals.append(np.sqrt(np.mean(data**2)))
#     residual_vals.append(data)
#     frame_num.append(i)   # 0 corresponds to 01124972
#
# table = pd.DataFrame({'Frame': frame_num, 'Mean': means, 'RMS': rms_vals, 'Residual': residual_vals})
# table.to_csv(r'D:\NLC\C1\frame_statistics.csv', index=False)

first_frame = 1124972

frame_num = []
means = []
rms_vals = []
median_vals = []
std_vals = []
for i in range(res.shape[0]):
    data = res[i]
    means.append(np.mean(data))
    rms_vals.append(np.sqrt(np.mean(data**2)))
    median_vals.append(np.median(data))
    std_vals.append(np.std(data))
    frame_num.append(first_frame + i)

table = pd.DataFrame({'Frame': frame_num, 'Mean': means, 'RMS': rms_vals, 'Median Res Val': median_vals, 'StdDev': std_vals})
table.to_csv(r'D:\NLC\C1\frame_statistics.csv', index=False)

std = np.nanstd(res)
for i in range(res.shape[0]):
    plt.figure()
    residuals_frame = res[i]
    bins = np.arange(-2*std, 2*std, std/20)
    hist = np.histogram(residuals_frame[np.isfinite(residuals_frame)], bins=bins)
    plt.bar(hist[1][:-1], hist[0], color='blue')
    plt.title(f'Histogram of Residuals for Frame {i + 1}')
    plt.xlabel('Residual Value')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.show()

