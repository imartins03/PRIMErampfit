from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

#definition of paths
super_bias_path = 'IRRC_calfiles\\super_biasC1.fits.ramp.20231012'
calFile = r'IRRC_calfiles\irrc_weights_C1.h5'
maskFile_path = r'IRRC_calfiles\C1_bad_ref_pix_mask.fits'
y_cube_path = r'D:\NLC\C1\y_cube_25.fits'
fit_cube_path = r'D:\NLC\C1\fit_cube.fits'
fit_coeff_path = r'D:\NLC\C1\fit_coeff.fits'
residuals_cube_path = r'D:\NLC\C1\residuals.fits'

#getting data from existing files
super_bias = fits.getdata(super_bias_path)  # Load super bias data
maskFile = fits.getdata(maskFile_path)  # Load mask file data
mask = maskFile > 0  # Create mask for bad pixels
supercpy = super_bias.copy()  # Create a copy of the super bias
supercpy[mask[:, :4096]] = 0  # Apply mask to the super bias copy

def evaluate_poly_array(coeffs, a_array, poly_type='power'):
    # Function to evaluate polynomial arrays
    output_arrays = []
    for a in a_array:  # Loop over input values
        print(a)
        if poly_type == 'power':  # Only 'power' type is considered here
            output_array = np.zeros(coeffs.shape[1])  # Initialize output array
            for n, coeff in enumerate(coeffs):  # Loop over coefficients
                output_array += coeff * (a ** n)  # Calculate polynomial value
            output_arrays.append(output_array)  # Append result to list
    return np.asarray(output_arrays)  # Convert list to numpy array

def generate_fit_cube(frame_num, degrees, saturation=50000, n_frames=None):
    y_cube = fits.getdata(y_cube_path)  # Load y_cube data
    print(np.ndim(y_cube))
    x = y_cube.shape[0]  # x is the dimension of the data cube (number of frames)

    y = y_cube.reshape(x, 4088, 4088)
    if n_frames is not None:
        y = y[:n_frames, :, :]  # Use specified number of frames
        x = n_frames

    print(y.shape)
    y = y_cube.reshape(x, -1)  # Reshape y_cube for fitting

    #saturation_mask = y > saturation
    # sat_pix = (saturation_mask).astype(float)


    print(y.shape)
    z = np.arange(len(y))  # Generate array for fitting
    print(z)

    coefficients, _ = np.polyfit(z, y, degrees, cov=True)  # Fit polynomial
    #coefficients[:, sat_pix.sum(axis=0) > 0] = np.nan

    print(coefficients.shape)
    fit_coeff = coefficients.reshape(degrees + 1, 4088, 4088)  # Reshape coefficients
    fits.writeto(fit_coeff_path, fit_coeff, overwrite=True)  # Save coefficients
    print(fit_coeff.shape)
    fit_cube = evaluate_poly_array(np.flip(coefficients, axis=0), z)  # Evaluate polynomial array
    print(y.shape[1])
    print(fit_cube.shape)

    fit_cube = fit_cube.reshape(len(z), 4088, 4088)  # Reshape fit cube
    fits.writeto(fit_cube_path, fit_cube, overwrite=True)  # Save fit cube
    return fit_cube

saturation = 50000 #currently not used
degrees = 6
generate_fit_cube(np.linspace(1, 25, 25), degrees, saturation)  # Generate fit cube

# def calculate_residuals():
#     y_cube = fits.getdata(y_cube_path)  # Load y_cube data
#     y_cube = y_cube[:, 0, :, :]  # Select the first element along the second axis
#     fit_cube = fits.getdata(fit_cube_path)  # Load fit cube data
#     residuals_cube = y_cube - fit_cube  # Calculate residuals
#     fits.writeto(residuals_cube_path, residuals_cube, overwrite=True)  # Save residuals cube
#     return residuals_cube
#
# res = calculate_residuals()  # Calculate residuals
#
# initial_frame_label = 1124972
#
# #creating lists to append values for the table
# frame_num = []
# means = []
# rms_vals = []
# median_vals = []
# std_vals = []
# for i in range(res.shape[0]):
#     data = res[i]
#     means.append(np.mean(data))  # Calculate mean
#     rms_vals.append(np.sqrt(np.mean(data ** 2)))  # Calculate RMS of residuals
#     median_vals.append(np.median(data))  # Calculate median
#     std_vals.append(np.std(data))  # Calculate std of residuals
#     frame_num.append(initial_frame_label + i)  # Adjusted frame numbering
#
# table = pd.DataFrame({'Mean': means, 'RMS': rms_vals, 'Median': median_vals, 'StdDev': std_vals})
# table.to_csv(r'D:\NLC\C1\frame_statistics.csv', index=False)  # Save statistics to CSV
#
# std = np.nanstd(res)
# for i in range(res.shape[0]):
#     residuals_frame = res[i]
#     bins = np.arange(-2 * std, 2 * std, std / 20)  # Generate histogram bins
#     hist = np.histogram(residuals_frame[np.isfinite(residuals_frame)], bins=bins)
#     plt.bar(hist[1][:-1], hist[0], color='blue')  # Plot histogram
#     plt.title(f'Histogram of Residuals for Frame {i + 1}')
#     plt.xlabel('Residual Value')
#     plt.ylabel('Frequency')
#     plt.grid(True)
#     plt.clf()
#     plt.savefig(f'D:\\NLC\\C1\\hist_{i}.png')  # Save histogram plot
