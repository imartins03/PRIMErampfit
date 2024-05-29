import numpy as np
from astropy.io import fits

# Constants
CUTOFF_GAIN = 5
ERR_THRESHOLD = 1125
FULL_WELL = 65535  # full-well (saturation) cutoff
DETECTOR_GAIN = 2
NUM_PIXELS = 524288  # Num Dat 128*4096 for an output

image1 = r'D:\NLC\C1\01124972C1_ircc.fits'
image2 = r'D:\NLC\C1\01124973C1_ircc.fits'
image3 = r'D:\NLC\C1\01124974C1_ircc.fits'
image4 = r'D:\NLC\C1\01124975C1_ircc.fits'

filenames = [image1, image2,image3,image4]
num_samples = len(filenames)

def read_fits_data(filenames):
    """Read and flatten data from a list of FITS files."""
    data_samples = []
    for f in filenames:
        with fits.open(f) as hdul:
            data_samples.append(hdul[0].data.reshape(-1))
    return np.array(data_samples)

def integrate(num_samples, sample_data, signal_out, variance_out):
    """Integrate sample data to compute signal and variance estimates."""
    num_pixels = sample_data.shape[1]

    # Initialize memory
    m = {
        'weight': np.zeros(num_pixels),
        'cumulative_signal': np.zeros(num_pixels),
        'total': np.zeros(num_pixels),
        'recent_sample': np.zeros(num_pixels),
        'sum_diff': np.zeros(num_pixels),
        'good_samples': np.zeros(num_pixels),
        'good_intervals': np.zeros(num_pixels),
    }

    # Process each sample
    for i in range(num_samples):
        next_sample = sample_data[i]
        diff = next_sample - m['recent_sample']
        avg_signal = m['cumulative_signal'] / np.maximum(m['weight'], 1)  # Avoid division by zero

        # Calculate error
        err = diff - avg_signal

        # Update memory based on conditions
        condition = sample_data[i] < FULL_WELL
        good_samples_condition = err * err < CUTOFF_GAIN * avg_signal + ERR_THRESHOLD

        m['weight'] += np.where(m['weight'] == 0, 1, 0)  # Set weight to 1 where it's 0
        m['cumulative_signal'] += np.where(good_samples_condition, m['good_samples']*next_sample-m['total'], 0)
        m['recent_sample'] = np.where(condition, next_sample, m['recent_sample'])
        m['total'] += np.where(condition, next_sample, 0)
        m['good_samples'] += np.where(good_samples_condition, 1, 0)
        m['good_intervals'] += np.where(condition, 1, 0)
        m['sum_diff'] += np.where(condition, diff, 0)

    # Calculate signal and variance estimates
    weight = m['weight']
    signal_out.append(m['cumulative_signal'] * (4 * ERR_THRESHOLD + m['sum_diff']) /
                      (m['good_intervals'] * m['cumulative_signal'] + 4 * weight * ERR_THRESHOLD))
    variance_out.append(ERR_THRESHOLD / (weight * DETECTOR_GAIN) +
                        signal_out[-1] / (m['good_intervals'] * DETECTOR_GAIN))

def main(input_files, output_file):
    """Main function to read data, integrate it, and save the result."""
    num_samples = len(input_files)
    sample_data = read_fits_data(input_files)
    signal_out = []
    variance_out = []

    integrate(num_samples, sample_data, signal_out, variance_out)

    with fits.open(input_files[0]) as hdul:
        hdul[0].data = np.array(signal_out).reshape(hdul[0].data.shape)
        output_file = 'output_ramp.fits'
        hdul.writeto(output_file, overwrite=True)

# how i would use it
input_files = [image1, image2]
output_file = 'output.ramp.fits'
main(input_files, output_file)