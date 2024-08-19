import numpy as np
import matplotlib.pyplot as plt

# Define the polynomial evaluation function for Legendre polynomials
def evaluate_legendre_poly(coeffs, x):
    """Evaluate Legendre polynomial with given coefficients."""
    return np.polynomial.legendre.legval(x, coeffs)

# Define specific polynomial coefficients (e.g., for a 10th-degree polynomial)
specific_coeffs = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])  # Coefficients in increasing order

# Parameters
degree = len(specific_coeffs) - 1
num_points = 100

# Generate synthetic data
x = np.linspace(0, 10, num_points)  # Change x to a time-like variable ranging from 0 to 10
y_true = evaluate_legendre_poly(np.flip(specific_coeffs), np.linspace(-1, 1, num_points))  # Generate true data

# Fit polynomial using Legendre basis
fit_coeffs = np.polynomial.legendre.legfit(np.linspace(-1, 1, num_points), y_true, degree)

# Evaluate the fitted polynomial using `legval`
y_fit = evaluate_legendre_poly(fit_coeffs, np.interp(x, [0, 10], [-1, 1]))  # Interpolate x to fit the Legendre basis range

# Compute residuals
residuals = y_true - y_fit

# Plot results
plt.figure(figsize=(12, 8))

# Plot true polynomial
plt.subplot(3, 1, 1)
plt.plot(x, y_true, 'r', label='True Polynomial')
plt.title('True Polynomial')
plt.legend()

# Plot fitted polynomial
plt.subplot(3, 1, 2)
plt.plot(x, y_fit, 'g', label='Fitted Polynomial')
plt.title('Fitted Polynomial')
plt.legend()

# Plot residuals
plt.subplot(3, 1, 3)
plt.scatter(x, residuals, label='Residuals')
plt.axhline(0, color='black', linestyle='--')
plt.title('Residuals')
plt.legend()

plt.tight_layout()
plt.show()
