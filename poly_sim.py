import numpy as np
import matplotlib.pyplot as plt

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

# Define specific polynomial coefficients (e.g., for a 10th-degree polynomial)
specific_coeffs = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11])  # Coefficients in increasing order

# Parameters
degree = len(specific_coeffs) - 1
num_points = 100

# Generate synthetic data
x = np.linspace(0, 10, num_points)
y_true = evaluate_poly_array(np.flip(specific_coeffs), x)  # Generate true data using specific polynomial

# Fit polynomial
fit_coeffs = np.polyfit(x, y_true, degree)

# Evaluate the fitted polynomial using your function
y_fit = evaluate_poly_array(np.flip(fit_coeffs), x)

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
