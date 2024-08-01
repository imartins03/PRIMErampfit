from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Paths
super_bias_path = 'IRRC_calfiles\\super_biasC1.fits.ramp.20231012'
maskFile_path = r'IRRC_calfiles\C1_bad_ref_pix_mask.fits'
y_cube_path = r'D:\NLC\C1\y_cube_100.fits'
fit_cube_base_path = r'D:\NLC\C1\fit_cube_'
fit_coeff_path = r'D:\NLC\C1\fit_coeff_'

# Load data
super_bias = fits.getdata(super_bias_path)  # Load super bias data
maskFile = fits.getdata(maskFile_path)  # Load mask file data
mask = maskFile > 0  # Create mask for bad pixels
supercpy = super_bias.copy()  # Create a copy of the super bias
supercpy[mask[:, :4096]] = 0  # Apply mask to the super bias copy
n_frames = 100

def evaluate_poly_array(coeffs, a_array, poly_type='power'):
    output_arrays = []
    for a in a_array:
        if poly_type == 'power':
            output_array = np.zeros(coeffs.shape[1])
            for n, coeff in enumerate(coeffs):
                output_array += coeff * (a ** n)
            output_arrays.append(output_array)
    return np.asarray(output_arrays)

def generate_fit_cube(center, degrees, saturation=50000):
    half_size = 128  # Since size=256, half_size is 128
    y_start = center[0] - half_size
    y_end = center[0] + half_size
    x_start = center[1] - half_size
    x_end = center[1] + half_size
    region = (y_start, y_end, x_start, x_end)

    y_cube = fits.getdata(y_cube_path)[1:n_frames]

    x, _, y, z = y_cube.shape
    y = y_cube.reshape(x, 4088, 4088)
    y = y[:, region[0]:region[1], region[2]:region[3]]
    y = y.reshape(x, -1)

    print('THE LENGTH OF Y IS', len(y))
    time = np.arange(len(y), dtype=np.double)
    print('y shape:', np.shape(y))

    # Perform polynomial fitting
    coefficients, _ = np.polyfit(time, y, degrees, cov=True)
    fit_coeff = coefficients.reshape(degrees + 1, region[1] - region[0], region[3] - region[2])

    # Debug the center value and filename strings
    print(f"Center coordinates: {center}")
    center_str = f"{center[1]}_{center[0]}"
    print(f"Center string: {center_str}")

    coeff_filename = f"{fit_coeff_path}center_{center_str}_{degrees}deg_noframe1.fits"
    fit_cube_filename = f"{fit_cube_base_path}center_{center_str}_{degrees}deg_noframe1.fits"

    # Debug the filenames
    print(f"Coefficient file: {coeff_filename}")
    print(f"Fit cube file: {fit_cube_filename}")

    # Save the FITS files
    fits.writeto(coeff_filename, fit_coeff, overwrite=True)

    fit_cube = evaluate_poly_array(np.flip(coefficients, axis=0), time)
    fit_cube = fit_cube.reshape(len(time), region[1] - region[0], region[3] - region[2])
    fits.writeto(fit_cube_filename, fit_cube, overwrite=True)

    return fit_cube

def process_superpixels(centers, size=256, degrees_range=(1, 10)):
    half_size = size // 2

    # Define regions
    regions = []
    for center in centers:
        y_start = center[0] - half_size
        y_end = center[0] + half_size
        x_start = center[1] - half_size
        x_end = center[1] + half_size
        regions.append((y_start, y_end, x_start, x_end))

    fit_results = {}
    for degree in range(degrees_range[0], degrees_range[1] + 1):
        for i, region in enumerate(regions):
            print(f"Processing region: {region} with degree: {degree}")  # Debugging statement
            fit_cube = generate_fit_cube(centers[i], degree)  # Ensure correct center and degree are passed
            fit_results[f'Region_{i}_Degree_{degree}'] = fit_cube
            print('fit_cube', fit_cube)
            print(np.shape(fit_cube))

    return fit_results

#this is just for my visualization purposes
def plot_superpixels_centers(image_size, centers, size=256):
    fig, ax = plt.subplots()

    # Create an empty image
    ax.imshow(np.zeros(image_size), cmap='gray', origin='lower')

    # Plot centers as dots and rectangles
    half_size = size // 2
    for center in centers:
        # Plot center as a red dot
        ax.plot(center[1], center[0], 'ro', markersize=10)  # Red dots for centers

        # Annotate the center with its coordinates
        ax.text(center[1] + half_size + 10, center[0] + half_size + 10,
                f'({center[1]}, {center[0]})', color='white', fontsize=8,
                verticalalignment='bottom', horizontalalignment='left')

        # Define the rectangle
        rect = patches.Rectangle(
            (center[1] - half_size, center[0] - half_size),
            size, size,
            linewidth=1, edgecolor='r', facecolor='none'
        )
        # Add the rectangle to the plot
        ax.add_patch(rect)

    ax.set_title('Superpixel Centers')
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.grid(True)
    plt.show()

centers = [(2048, 2048), (3072, 2048), (500, 2048)]  # Centers of the superpixels
degrees_range = (1, 10)  # Range of degrees

# Process superpixels
fit_results = process_superpixels(centers, size=256, degrees_range=degrees_range)

# Plot the centers of superpixels
plot_superpixels_centers(image_size=(4088, 4088), centers=centers, size=256)
