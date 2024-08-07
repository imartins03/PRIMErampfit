from astropy.io import fits
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# Definition of paths
super_bias_path = 'IRRC_calfiles\\super_biasC1.fits.ramp.20231012'
calFile = r'IRRC_calfiles\irrc_weights_C1.h5'
maskFile_path = r'IRRC_calfiles\C1_bad_ref_pix_mask.fits'
# y_cube_path = r'D:\NLC\C1\y_cube_100.fits'
# stat_table_path = r'D:\NLC\C1\avg_signal_y_cube.csv'
#
# y_cube = fits.getdata(y_cube_path)
#
# def compute_statistics(y_cube, initial_frame_label):
#     means = []
#     frame_num = []
#
#     for i in range(y_cube.shape[0]):
#         data = y_cube[i]
#         means.append(np.mean(data))  # Calculate mean
#         frame_num.append(initial_frame_label + i)  # Adjusted frame numbering
#
#     table = pd.DataFrame({
#         'FrameNumber': frame_num,
#         'Average Signal': means,
#     })
#
#     return table
#
# initial_frame_label = 1124972
# statistics_df = compute_statistics(y_cube, initial_frame_label)  # Compute statistics
#
# # Save the DataFrame to a CSV file
# statistics_df.to_csv(stat_table_path, index=False)
#
# print(f'Statistics table saved to {stat_table_path}')

#%%

# from astropy.io import fits
# import numpy as np
# import pandas as pd
#
# stat_table_path_raw = r'D:\NLC\C1\avg_signal_raw_frames.csv'
# file_format = r'D:\NLC\C1\{0:08d}C1.fits.fz'
# frame_start = 1124972
# frame_count = 100
# file_list = [file_format.format(frame_start+i) for i in range(frame_count)] #method from before did not work, was doing 40000 frames
#
# def compute_statistics(file_list, initial_frame_label):
#     means = []
#     frame_num = []
#
#     # Iterate over each file
#     for file_path in file_list:
#         # Read data from the FITS file
#         data = fits.getdata(file_path)
#
#
#         # Check the number of frames in the file and ensure we're processing only up to the limit
#         num_frames = data.shape[0]
#
#         # Process frames up to the frame_count
#         for i in range(min(num_frames, frame_count)):
#             frame_data = data[i]
#             means.append(np.mean(frame_data))  # Calculate mean
#             frame_num.append(initial_frame_label)
#             initial_frame_label += 1  # Increment frame label for each new frame
#
#             # Stop processing if we have reached the desired frame_count
#             if len(means) >= frame_count:
#                 break
#
#         if len(means) >= frame_count:
#             break
#
#     # Create a DataFrame to store the results
#     table = pd.DataFrame({
#         'FrameNumber': frame_num,
#         'Average Signal': means,
#     })
#
#     return table
#
#
# # Initial frame label for the processing
# initial_frame_label = frame_start
# statistics_df = compute_statistics(file_list, initial_frame_label)
#
# # Save the DataFrame to a CSV file
# statistics_df.to_csv(stat_table_path_raw, index=False)
#
# print(f'Statistics table saved to {stat_table_path_raw}')

#%%
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

# File path
y_cube_path = r'D:\NLC\C1\y_cube_500.fits'

def compute_statistics(y_cube_path, initial_frame_label):
    means = []
    frame_num = []
#
    # Read data from the FITS file
    y_cube = fits.getdata(y_cube_path)

    # Iterate over each frame in the cube
    for i in range(y_cube.shape[0]):
        frame_data = y_cube[i]
        means.append(np.mean(frame_data))  # Calculate mean ADU value for the frame

    return means

# Initial frame label for the processing
initial_frame_label = 1124972

# Compute statistics
average_adus = compute_statistics(y_cube_path, initial_frame_label)

# Plot the average ADU value as a function of frame number
plt.figure(figsize=(10, 6))
plt.plot(np.arange(len(average_adus)), average_adus, marker='o', linestyle='-')
plt.xlabel('Frame Number')
plt.ylabel('Average ADU Value')
plt.title('Average ADU Value per Frame')
plt.grid(True)
plt.savefig(r'D:\NLC\C1\average_adu_per_frame.png')  # Save the plot as an image file
plt.show()

saturation_val = 50000
y_cube = fits.getdata(y_cube_path)
print(y_cube.shape)
y = y_cube[:,0,:,:]

def percentage_saturation_val():
    percentages = []

    for i in range(y.shape[0]):
        frame_data = y[i]

        pix_over_saturation = np.sum(frame_data > saturation_val)

        # Calculate the total number of pixels
        total_pixels = 4088*4088

        # Calculate the percentage of saturated pixels
        percentage_pix_over_saturation = (pix_over_saturation/total_pixels)*100
        percentages.append(percentage_pix_over_saturation)
        print(f"Frame {i}: {percentage_pix_over_saturation:.2f}%")

    return percentages

percentages = percentage_saturation_val()



