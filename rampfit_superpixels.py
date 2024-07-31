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

    y_cube = fits.getdata(y_cube_path)[1:n_frames]   #I did this slicing so that it cuts out the first frame

    x,_, y, z = y_cube.shape
    # Extract the region
    y = y_cube.reshape(x, 4088, 4088)
    y = y[:, region[0]:region[1], region[2]:region[3]]  #this defines where on the image the superpixel is

    # Reshape for fitting
    y = y.reshape(x, -1)  # Flatten the spatial dimensions, this makes the shape of y (4,65536) which is (4, 256*256)

    print('THE LENGTH OF Y IS', len(y))
    time = np.arange(len(y), dtype=np.double) #len(y) is still just the number of frames
    print('y shape:',np.shape(y))

    # Perform polynomial fitting
    coefficients, _ = np.polyfit(time, y, degrees, cov=True)
    fit_coeff = coefficients.reshape(degrees + 1, region[1] - region[0], region[3] - region[2]) #region[1]-region[0] in this case will be 256 because thats the size of the superpixel (degrees+1,256,256)
    center_str = f"{center[0]}_{center[1]}"
    fits.writeto(fit_coeff_path + f'center_{center_str}_{degrees}_noframe1.fits', fit_coeff, overwrite=True)

    fit_cube = evaluate_poly_array(np.flip(coefficients, axis=0), time)
    fit_cube = fit_cube.reshape(len(time), region[1] - region[0], region[3] - region[2])
    fits.writeto(fit_cube_base_path + f'center_{center_str}_{degrees}_noframe1.fits', fit_cube, overwrite=True)

    return fit_cube

def process_superpixels(centers, size=256, degrees=1):
    fit_results = {}
    for center in centers:
        fit_cube = generate_fit_cube(center, degrees)
        center_str = f"{center[0]}_{center[1]}"
        fit_results[f'Center_{center_str}'] = fit_cube
        print(f'fit_cube for center {center}:', fit_cube)

    return fit_results


#I just did this to give me a visual representation of where the centers of the three superpixels are
def plot_superpixels_centers(image_size, centers, size=256):
    fig, ax = plt.subplots()
    ax.imshow(np.zeros(image_size), cmap='gray', origin='lower')

    half_size = size // 2

    # Plot centers as dots and squares
    for center in centers:
        # Plot the center point
        ax.plot(center[1], center[0], 'ro', markersize=10)  # Red dots for centers

        # Define the square's position
        rect = patches.Rectangle(
            (center[1] - half_size, center[0] - half_size),  # Bottom-left corner
            size,  # Width
            size,  # Height
            linewidth=1,  # Line width
            edgecolor='r',  # Edge color
            facecolor='none'  # No fill color
        )
        # Add the rectangle to the plot
        ax.add_patch(rect)

    ax.set_title('Superpixel Centers')
    ax.set_xlabel('X Coordinate')
    ax.set_ylabel('Y Coordinate')
    ax.grid(True)
    plt.show()

centers = [(2048, 2048), (3072,2048), (1024,2048)]  # Centers of the superpixels
degrees = 6

# Process superpixels
fit_results = process_superpixels(centers, size=256, degrees=degrees)

# Plot the centers of superpixels
plot_superpixels_centers(image_size=(4088, 4088), centers=centers, size=256)
