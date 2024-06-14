
import os
import numpy as np
from astropy.io import fits
import pandas as pd


fuv = pd.read_csv("every3_fuv.csv")
eclipse_list = fuv['eclipse'][0:]

cumulative_image = None
eclipse_counter = 0
file = 0
for e in eclipse_list:
    if eclipse_counter <= 2500:
        i = str(e).zfill(5)
        # for mounted backplanetest bucket
        file_path = f"/mnt/s3/e{i}-fd--NF/e{i}-fd-b00-f0100-t00000-g_dose.fits.gz"

        if os.path.isfile(file_path):
            try:
                with fits.open(file_path) as hdul:
                    data = hdul[0].data
                    if data is None:
                        print("no data")
                    if cumulative_image is None:
                        cumulative_image = np.zeros_like(data, dtype=np.float64)
                    if data is not None:
                        print(f"adding eclipse {e}, spot {eclipse_counter}")
                        data[data > 2] = 0
                        cumulative_image += data.astype(np.float64)
                        eclipse_counter = eclipse_counter + 1
            except Exception as ex:
                print(f"probably a corrupt fits file, {ex}")
        else:
            print(f"dose map ({file_path}) does not exist")

hdu = fits.PrimaryHDU(cumulative_image)
hdul = fits.HDUList([hdu])
combo_filename = "/mnt/s3/fuv_2500stack_filtered.fits"
hdul.writeto(combo_filename, overwrite=True)

print(f"Stacked FITS file '{combo_filename}' has been saved.")
