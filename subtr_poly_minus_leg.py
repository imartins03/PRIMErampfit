import numpy as np
from astropy.io import fits

coeff_deg1_239_poly = fits.getdata(r'F:\leftover_C1_dif_degrees_test_rampfit\239_frames\unweighted\fit_coeff_poly_1deg_239frames_noframe1.fits')
coeff_deg1_239_leg = fits.getdata(r"F:\legfit\239_frames_unweighted\coefficients_leg_239frames_1deg_final_vst_unweighted.fits")

def coeff_res():
    res = coeff_deg1_239_poly - coeff_deg1_239_leg
    print(res)
    return res

coeff_res()