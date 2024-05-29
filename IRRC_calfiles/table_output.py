import os
import pandas as pd
from astropy.io import fits

# Directory containing the data files
data_directory = "C:\PycharmProjects\PRIMErampfit"

all_data = []

for filename in os.listdir(data_directory):
    if filename.endswith("ramp.fits"):  #??
        filepath = os.path.join(data_directory, filename)
        with fits.open(filepath) as hdul:
            header = hdul[0].header
            # Extracting data from the header and formatting it

            file_data = {
                "DATEBEG": header.get('DATE-BEG'),
                "filename": filename,
                "OBJNAME": header.get('OBJNAME'),
                "OBJTYPE": header.get('OBJTYPE'),
                "OBSERVER": header.get('OBSERVER'),
                "CHIP": header.get('CHIP'),
                "FILTER1": header.get('FILTER1'),
                "FILTER2": header.get('FILTER2'),
                "EXPTIME": header.get('EXPTIME'),
                "NFRAMES": header.get('NFRAMES'),
                "DITHTYP": header.get('DITHTYP'),
                "DITHRAD": header.get('DITHRAD'),
                "DITHPH": header.get('DITHPH'),
                "DITH_REP": header.get('DITH_REP'),
                "AIRMASS": header.get('TSSECZ','N/A'),
                "FOCUS": header.get('FOCUS'),
                "NINT": header.get('NINT'),
                "INT": header.get('INT')
            }
            all_data.append(file_data)


df = pd.DataFrame(all_data)
output_csv = 'output_data.csv'
df.to_csv(output_csv, index=False)