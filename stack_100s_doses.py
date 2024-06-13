
import os
import numpy as np
from astropy.io import fits
import pandas as pd


fuv = pd.read_csv("every3_fuv.csv")
eclipse_list = fuv['eclipse'][0:]

cumulative_image = None
eclipse_counter = 0

for e in eclipse_list:
    while eclipse_counter <= 2000:
        i = str(e).zfill(5)
        # for mounted backplanetest bucket
        file_path = f"/mnt/s3/e{i}-fd--NF/e{i}-fd-b00-f0100-t00000-g_dose.fits.gz"

        if os.path.isfile(file_path):
            with fits.open(file_path) as hdul:
                data = hdul[0].data
                if data is None:
                    print("no data")
                if cumulative_image is None:
                    cumulative_image = np.zeros_like(data, dtype=np.float64)
                if data is not None:
                    cumulative_image += data.astype(np.float64)
                    eclipse_counter = eclipse_counter + 1
        else:
            print(f"dose map ({file_path}) does not exist")

hdu = fits.PrimaryHDU(cumulative_image)
hdul = fits.HDUList([hdu])
combo_filename = "/mnt/s3/fuv_2000stack_100seceach.fits"
hdul.writeto(combo_filename, overwrite=True)

print(f"Stacked FITS file '{combo_filename}' has been saved.")
