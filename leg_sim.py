import numpy as np
import matplotlib.pyplot as plt
import os

def evaluate_poly_array(coeffs, a_array, poly_type='power'):
    output_arrays = []
    for a in a_array:
        if poly_type == 'power':
            output_array = 0
            for n, coeff in enumerate(coeffs):
                output_array += coeff * (a ** n)
            output_arrays.append(output_array)
    return np.array(output_arrays)

# Define the polynomial evaluation function for Legendre polynomials
def evaluate_legendre_poly(coeffs, x):
    """Evaluate Legendre polynomial with given coefficients."""
    return np.transpose(np.polynomial.legendre.legval(x, coeffs))

def save_plot(x, y_true, y_fit, residuals, degree, output_dir):
    """Save plots of true polynomial, fitted polynomial, and residuals."""
    plt.figure(figsize=(12, 8))

    # Plot true polynomial
    plt.subplot(3, 1, 1)
    plt.plot(x, y_true, 'r', label='True Polynomial')
    plt.title(f'True Polynomial (Degree {degree})')
    plt.legend()

    # Plot fitted polynomial
    plt.subplot(3, 1, 2)
    plt.plot(x, y_fit, 'g', label='Fitted Polynomial')
    plt.title(f'Fitted Legendre Poly (Degree {degree})')
    plt.legend()

    # Plot residuals
    plt.subplot(3, 1, 3)
    plt.scatter(x, residuals, label='Residuals', color='blue')
    plt.axhline(0, color='black', linestyle='--')
    plt.title('Residuals')

    # Save figure
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'legendre_fit_degree_{degree}.png'))
    plt.show()
    plt.close()

def generate_fit_plots(degree, num_points=100):
    """Generate and save plots for a Legendre polynomial fit of a specified degree."""
    # Generate synthetic data
    x = np.linspace(0, 100, num_points+1)
    # x = np.linspace(-1, 1, num_points)
    specific_coeffs = np.arange(degree + 1)  # Use coefficients from 0 to degree

    y_true = evaluate_poly_array(np.flip(specific_coeffs, axis=0), x)  # Evaluate polynomial array

    # Fit polynomial using Legendre basis
    fit_coeffs = np.polynomial.legendre.legfit(x, y_true, degree)

    # Evaluate the fitted polynomial
    y_fit = evaluate_legendre_poly(fit_coeffs, x)

    # Compute residuals
    residuals = y_true - y_fit

    # Create output directory
    output_dir = f'F:'

    # Save plots
    save_plot(x, y_true, y_fit, residuals, degree, output_dir)

# Generate plots for Legendre polynomial of degree 10
generate_fit_plots(degree=10)
