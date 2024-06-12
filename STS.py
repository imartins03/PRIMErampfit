import calibnonlinCorrection as cal
from calib_nonlinCorrection import *
from Demos.BackupRead_BackupWrite import outfile
import math
import sys
import os
import glob
from astropy.io import fits
import numpy as np
from irrc.apply_correction import _load_coeffs_from_file, apply_in_memory
from irrc.util import destripe
import matplotlib.pyplot as plt
import pandas as pd
import scipy
from scipy.optimize import curve_fit
from astropy.modeling.polynomial import Polynomial1D, Legendre1D
from astropy.modeling import fitting
import string

#seems essentially the same to get_ramp_cds
def get_ramp(file,pixX,pixY):
    '''PixX is first coordinate of pixel,pixY is second coordinate of pixel
    RETUNRNS array with pixel coordinates
    '''