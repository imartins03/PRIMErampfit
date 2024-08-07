from astropy.io import fits
import numpy as np

# Template path for fit coefficient FITS file
fit_coeff_path_template = r'F:\leftover_C1_dif_degrees_test_rampfit\239_frames\fit_coeff_poly_{degree}deg_239frames_noframe1.fits'


def load_coefficients(degree):
    fit_coeff_path = fit_coeff_path_template.format(degree=degree)
    fit_coeff = fits.getdata(fit_coeff_path)
    return fit_coeff


def output_polynomial(coefficients):
    degree = len(coefficients) - 1
    poly_terms = []
    for n in range(degree + 1):
        coeff = coefficients[n]
        if n == 0:
            term = f"{coeff:.6e}"
        else:
            term = f"{coeff:.6e} * x^{n}"
        poly_terms.append(term)

    polynomial = " + ".join(poly_terms)
    print(f"Polynomial of degree {degree}:")
    print(polynomial)


if __name__ == "__main__":
    # Example usage
    degree = 3  # Example degree (you can choose any degree you want)
    coefficients = load_coefficients(degree)
    output_polynomial(coefficients)
