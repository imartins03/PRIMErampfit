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
# y_cube_path = r'D:\NLC\C1\y_cube_100.fits'
stat_table_path_raw = r'D:\NLC\C1\avg_signal_raw_frames.csv'

file_format = r'D:\NLC\C1\{0:08d}C1.fits.fz'
frame_range = (1124972, 1124972 + 100)
file_list = [file_format.format(n) for n in range(*frame_range)]

def compute_statistics(file_list, initial_frame_label):
    means = []
    frame_num = []

    for file_path in file_list:
        data = fits.getdata(file_path)
        for i in range(data.shape[0]):
            frame_data = data[i]
            means.append(np.mean(frame_data))  # Calculate mean
            frame_num.append(initial_frame_label)
            initial_frame_label += 1  # Increment frame label for each new frame

    table = pd.DataFrame({
        'FrameNumber': frame_num,
        'Average Signal': means,
    })

    return table

initial_frame_label = 1124972
statistics_df = compute_statistics(file_list, initial_frame_label)

statistics_df.to_csv(stat_table_path_raw, index=False)