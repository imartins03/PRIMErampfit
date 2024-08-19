import numpy as np
import matplotlib.pyplot as plt

# Define your polynomial evaluation function
def evaluate_poly_array(coeffs, a_array, poly_type='power'):
    """Evaluate a polynomial with given coefficients at the points in a_array."""
    output_arrays = []
    for a in a_array:
        if poly_type == 'power':
            output_array = sum(coeff * (a ** n) for n, coeff in enumerate(coeffs))
            output_arrays.append(output_array)
    return np.array(output_arrays)

def generate_fit_plots(degree, num_points=100):
    """Generate residuals for a polynomial fit of a specified degree."""
    # Generate synthetic data
    x = np.linspace(1, num_points, num_points)
    specific_coeffs = np.arange(degree + 1)  # Coefficients from 0 to degree
    y_true = evaluate_poly_array(np.flip(specific_coeffs, axis=0), x)

    # Fit polynomial using np.polyfit
    fit_coeffs = np.polyfit(x, y_true, degree)

    # Evaluate the fitted polynomial
    y_fit = evaluate_poly_array(np.flip(fit_coeffs, axis=0), x)

    # Compute residuals
    residuals = y_true - y_fit

    # Return residuals
    return residuals

# Generate and plot median residuals for polynomial degrees 1 through 10
degrees = range(1, 11)
median_residuals = []

for degree in degrees:
    residuals = generate_fit_plots(degree)
    median_residuals.append(np.median(residuals))

plt.figure(figsize=(8, 6))
plt.plot(degrees, median_residuals, 'o-', label='Median Residual')
plt.xlabel('Polynomial Degree')
plt.ylabel('Median Residual')
plt.title('Median Residual vs Polynomial Degree')
plt.legend()
plt.grid(True)
plt.savefig(r'F:\median_residuals_poly.png')  # Use raw string for Windows path
plt.show()
