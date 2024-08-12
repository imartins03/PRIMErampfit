from astropy.io import fits
import numpy as np
import os

# Definition of paths
super_bias_path = 'IRRC_calfiles/super_biasC1.fits.ramp.20231012'
calFile = r'IRRC_calfiles/irrc_weights_C1.h5'
maskFile_path = r'IRRC_calfiles/C1_bad_ref_pix_mask.fits'
y_cube_path = r'D:\NLC\C1\y_cube_500.fits'


def evaluate_legendre_poly(coeffs, x):
    """Evaluate Legendre polynomial with given coefficients."""
    return np.polynomial.legendre.legval(x, coeffs)

def save_coefficients(coefficients, degree, row, output_dir):
    """Save Legendre polynomial coefficients to a .npy file."""
    coeffs_path = os.path.join(output_dir, f'coefficients_degree{degree}_row{row}.npy')
    np.save(coeffs_path, coefficients)


def save_fit_cube_row(fit_cube_row, degree, row, output_dir):
    """Save fit cube row to a FITS file."""
    fit_cube_row_path = os.path.join(output_dir, f'fit_cube_row{row}_degree{degree}.fits')
    fits.writeto(fit_cube_row_path, fit_cube_row, overwrite=True)


def generate_fit_cube(degree, saturation=50000, n_frames=20, num_rows=4):
    """Generate and save fit cubes for a given Legendre polynomial degree."""
    # Load y_cube data
    y_cube = fits.getdata(y_cube_path)

    # Ensure we have enough frames and reshape if necessary
    y_cube = y_cube[1:n_frames]  # Use specified number of frames

    x = y_cube.shape[0]
    y = y_cube.reshape(x, 4088, 4088)

    # Print shapes to verify
    print(f"Data shape after slicing and reshaping: {y.shape}")

    row_size = y.shape[1] // num_rows  # Calculate size of each row
    fit_cube = np.zeros((y_cube.shape[0], 4088, 4088), dtype=np.double)  # Initialize fit cube

    output_dir = f'F:/legfit'
    os.makedirs(output_dir, exist_ok=True)  # Create output directory if it doesn't exist

    # Process each row separately
    for i in range(num_rows):
        row_start = i * row_size
        row_end = (i + 1) * row_size if i < num_rows - 1 else y.shape[1]
        y_row = y[:, row_start:row_end, :]
        y_row = y_row.reshape(y_cube.shape[0], -1)  # Reshape for fitting

        time = np.linspace(-1, 1, y_cube.shape[0], dtype=np.double)  # Legendre time

        # Fit Legendre polynomial
        coefficients = np.polynomial.legendre.legfit(time, y_row, degree)
        print(f'Fitting coefficients for row {i}:', coefficients)

        # Save coefficients
        save_coefficients(coefficients, degree, i, output_dir)

        # Generate and save fit cube for this row
        fit_cube_row = evaluate_legendre_poly(coefficients, time)
        fit_cube_row = fit_cube_row.reshape(y_cube.shape[0], row_end - row_start, 4088)  # Reshape fit cube for this row

        # Save individual row fit cube
        save_fit_cube_row(fit_cube_row, degree, i, output_dir)

        # Place fit_cube_row in the correct place in the final fit_cube
        fit_cube[:, row_start:row_end, :] = fit_cube_row
        print(f'Fit cube shape for row {i}:', fit_cube_row.shape)

    # Save final stacked fit cube
    fit_cube_path = os.path.join(output_dir, f'fit_cube_leg_{degree}deg_final_vst.fits')
    fits.writeto(fit_cube_path, fit_cube, overwrite=True)  # Save final fit cube
    print('Final fit cube shape:', fit_cube.shape)

# Loop through degrees from 1 to 10
for degree in range(1, 11):
    generate_fit_cube(degree)
