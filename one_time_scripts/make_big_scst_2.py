from gPhoton.io.mast import download_data
from pyarrow import parquet
import pyarrow as pa
import pandas as pd 
import numpy as np 
from astropy.io import fits 
import random
from astropy.table import Table 
import os 
from astropy.table import vstack 


# get schema 

print("getting schema")
DESIRABLE_FIELDS = ['pktime','flags','ptag','eclipse','hvnom_nuv','hvnom_fuv','ra_acs','dec_acs','roll_acs'] #'ra_rta','dec_rta','roll_rta','status_flag_rta','acs_mode','gyro_rate_rss','comment']

file_path = download_data(10298,'scst')
eclipse = Table.read(file_path).to_pandas()
all_fields = {}
tab = pa.Table.from_pandas(eclipse)
all_fields |= {
        name: tab.schema.field(name).type
        for name in tab.column_names
    }

filtered_fields = {k: v for k, v in all_fields.items() if k in DESIRABLE_FIELDS}
full_schema = pa.schema(filtered_fields)
print("got schema")


writer = parquet.ParquetWriter('extended_scst_7_16.parquet', full_schema)

# second pass: actually make the table
parquet_writer = None # define the parquet writer

eclipse_group, eclipse_count, eclipse_threshold = [], 0, 1000

print("loading eclipse list")
eclipses = np.loadtxt('eclipseList.csv')
#eclipse_list = random.choices(eclipses, k=20)

for eclipse in eclipses:
    
    # eclipse 
    file_path = download_data(eclipse,'scst')
    table_fits = Table.read(file_path)
    eclipse_df = table_fits.to_pandas()[DESIRABLE_FIELDS]
    
    # add number of legs for [0] header 
    fits_header = fits.open(file_path)
    eclipse_df['MPSNPOS'] = fits_header[0].header['MPSNPOS']
    fits_header.close()
    
    # table clean up 
    eclipse_df = eclipse_df.drop(columns=[c for c in eclipse_df.columns if c not in DESIRABLE_FIELDS])
    for k in filtered_fields.keys():
        if k not in eclipse_df.columns:
            eclipse_df[k] = np.nan
            
    tab = pa.Table.from_pandas(eclipse_df)
    
    # clean up 
    os.remove(file_path)
    del eclipse_df
    
    tab = tab.cast(full_schema)
    tab = tab.replace_schema_metadata({})
    eclipse_group.append(tab)
    if len(eclipse_group) == eclipse_threshold or eclipse_count+1 == len(eclipses):
        # add pa tables for 1000 eclipses to parquet, written as a group 
        writer.write_table(pa.concat_tables(eclipse_group))
        eclipse_group = []
    eclipse_count = eclipse_count + 1 
    
    print(f"on eclipse # {eclipse_count}")
    
writer.close()
