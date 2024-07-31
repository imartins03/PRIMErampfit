from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# Paths
super_bias_path = 'IRRC_calfiles\\super_biasC1.fits.ramp.20231012'
maskFile_path = r'IRRC_calfiles\C1_bad_ref_pix_mask.fits'
y_cube_path = r'D:\NLC\C1\y_cube_4.fits'
fit_cube_base_path = r'D:\NLC\C1\fit_cube_region_'
fit_coeff_base_path = r'D:\NLC\C1\fit_coeff_region_'

# Load data
super_bias = fits.getdata(super_bias_path)  # Load super bias data
maskFile = fits.getdata(maskFile_path)  # Load mask file data
mask = maskFile > 0  # Create mask for bad pixels
supercpy = super_bias.copy()  # Create a copy of the super bias
supercpy[mask[:, :4096]] = 0  # Apply mask to the super bias copy


def evaluate_poly_array(coeffs, a_array, poly_type='power'):
    output_arrays = []
    for a in a_array:
        if poly_type == 'power':
            output_array = np.zeros(coeffs.shape[1])
            for n, coeff in enumerate(coeffs):
                output_array += coeff * (a ** n)
            output_arrays.append(output_array)
    return np.asarray(output_arrays)


def generate_fit_cube(region, degrees, saturation=50000):
    y_cube = fits.getdata(y_cube_path)
    x,_, y, z = y_cube.shape
    # Extract the region
    y = y_cube.reshape(x, 4088, 4088)
    y = y[:, region[0]:region[1], region[2]:region[3]]

    # Reshape for fitting
    y = y.reshape(x, -1)  # Flatten the spatial dimensions

    time = np.arange(len(y), dtype=np.double)

    # Perform polynomial fitting
    coefficients, _ = np.polyfit(time, y, degrees, cov=True)
    fit_coeff = coefficients.reshape(degrees + 1, region[1] - region[0], region[3] - region[2])
    fits.writeto(fit_coeff_base_path + f'region_{region[0]}_{region[2]}.fits', fit_coeff, overwrite=True)

    fit_cube = evaluate_poly_array(np.flip(coefficients, axis=0), time)
    fit_cube = fit_cube.reshape(len(time), region[1] - region[0], region[3] - region[2])
    fits.writeto(fit_cube_base_path + f'region_{region[0]}_{region[2]}.fits', fit_cube, overwrite=True)
    return fit_cube


def process_superpixels(centers, size=256, degrees=1):
    half_size = size // 2

    # Define regions relative to the provided center coordinates
    regions = []
    for center in centers:
        y_start = center[0] - half_size
        y_end = center[0] + half_size
        x_start = center[1] - half_size
        x_end = center[1] + half_size
        regions.append((y_start, y_end, x_start, x_end))

    fit_results = {}
    for i, region in enumerate(regions):
        fit_cube = generate_fit_cube(region, degrees)
        fit_results[f'Region_{i}'] = fit_cube

    return fit_results

def plot_superpixels_centers(image_size, centers, size=256):
    plt.figure(figsize=(10, 10))
    plt.imshow(np.zeros(image_size), cmap='gray', origin='lower')

    # Plot centers as dots
    for center in centers:
        plt.plot(center[1], center[0], 'ro', markersize=10)  # Red dots for centers

    plt.title('Superpixel Centers')
    plt.xlabel('X Coordinate')
    plt.ylabel('Y Coordinate')
    plt.grid(True)
    plt.show()


centers = [(2048, 2048), (3072,2048), (1024,2048)]  # Centers of the superpixels
degrees = 1

# Process superpixels
fit_results = process_superpixels(centers, size=256, degrees=degrees)

# Plot the centers of superpixels
plot_superpixels_centers(image_size=(4088, 4088), centers=centers, size=256)
