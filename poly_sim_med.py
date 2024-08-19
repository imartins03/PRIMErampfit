import numpy as np
import matplotlib.pyplot as plt
import os

# Define your polynomial evaluation function
def evaluate_poly_array(coeffs, a_array, poly_type='power'):
    output_arrays = []
    for a in a_array:
        if poly_type == 'power':
            output_array = 0
            for n, coeff in enumerate(coeffs):
                output_array += coeff * (a ** n)
            output_arrays.append(output_array)
    return np.array(output_arrays)

# def save_plot(x, y_true, y_fit, residuals, degree, output_dir):
#     """Save plots of true polynomial, fitted polynomial, and residuals."""
#     plt.figure(figsize=(12, 8))
#
#     # Plot true polynomial
#     plt.subplot(3, 1, 1)
#     plt.plot(x, y_true, 'r', label='True Polynomial')
#     plt.title(f'True Polynomial (Degree {degree})')
#     plt.legend()
#
#     # Plot fitted polynomial
#     plt.subplot(3, 1, 2)
#     plt.plot(x, y_fit, 'g', label='Fitted Polynomial')
#     plt.title(f'Fitted Polynomial (Degree {degree})')
#     plt.legend()
#
#     # Plot residuals
#     plt.subplot(3, 1, 3)
#     plt.scatter(x, residuals, label='Residuals', color='blue')
#     plt.axhline(0, color='black', linestyle='--')
#     plt.title('Residuals')
#
#     # Save figure
#     plt.tight_layout()
#     plt.savefig(os.path.join(output_dir, f'poly_fit_degree_{degree}.png'))
#     plt.show()
#     plt.close()

def generate_fit_plots(degree, num_points=100):
    """Generate and save plots for a polynomial fit of a specified degree."""
    # Generate synthetic data
    x = np.linspace(1, 100, num_points)
    specific_coeffs = np.arange(degree + 1)  # Use coefficients from 0 to degree
    y_true = evaluate_poly_array(np.flip(specific_coeffs, axis=0), x)

    # Fit polynomial using np.polyfit
    fit_coeffs = np.polyfit(x, y_true, degree)

    # Evaluate the fitted polynomial
    y_fit = evaluate_poly_array(np.flip(fit_coeffs, axis=0), x)

    # Compute residuals
    residuals = y_true - y_fit

    # Save plots
    # save_plot(x, y_true, y_fit, residuals, degree, output_dir)
    return np.median(residuals)

# Generate and plot median residuals for polynomial degrees 1 through 10
degrees = range(1, 11)
median_residuals = [generate_fit_plots(degree) for degree in degrees]

plt.figure(figsize=(8, 6))
plt.plot(degrees, median_residuals, 'o-', label='Median Residual')
plt.xlabel('Polynomial Degree')
plt.ylabel('Median Residual')
plt.title('Median Residual vs Polynomial Degree')
plt.legend()
plt.grid(True)
plt.savefig('F:\median_residuals_poly.png')
plt.show()
