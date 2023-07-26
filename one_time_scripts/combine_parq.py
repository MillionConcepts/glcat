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


scst = parquet.read_table("extended_scst_7_16.parquet",columns=['eclipse','pktime','hvnom_nuv']).to_pandas()

scst = scst.rename(columns={"pktime": "time"})

parq = parquet.read_table("aspect.parquet",columns=['eclipse','time'])
aspect = parq.to_pandas()
aspect['refined'] = 1 
    
# merge to get slew frames
slew_frames = scst.merge(aspect, how='left', on=['time'])
slew_frames = slew_frames[slew_frames['hvnom_nuv'] == 1]

# use acs solns for ra / dec / roll for slew frames in asp soln
#slew_frames['ra_combo'] = slew_frames['ra'].fillna(slew_frames['ra_acs'])
#slew_frames['dec_combo'] = slew_frames['dec'].fillna(slew_frames['dec_acs'])
#slew_frames['roll_combo'] = slew_frames['roll'].fillna(slew_frames['roll_acs'])
#slew_frames = slew_frames.astype({
#                          'ra': 'float64',
#                          'dec': 'float64',
#                          'roll': 'float64'})

# need to add a column to indicate slew frames in the final aspect table (bc ra and dec are getting modded)
# save to parquet
slew_frames.to_parquet("aspect_with_slew.parquet", compression='snappy')
