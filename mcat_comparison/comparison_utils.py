import matplotlib.pyplot as plt
from pyarrow import parquet as pq
import pyarrow as pa
import pandas as pd
from hostess.aws.s3 import Bucket
import os
import pdr
import numpy as np
from gPhoton.pipeline import execute_pipeline

import warnings
# Suppress all UserWarnings
warnings.filterwarnings('ignore', category=UserWarning)

def get_mcat_file(eclipse,rootpath,
                  visit_metadata_fn='../glcat_v1_pipeline/visit_metadata.parquet',
                  mcat_manifest_fn = '../glcat_v1_pipeline/mcat_mast_list.csv',
                  bucket_name = 'uraniborg-sieve-7738937'):
    # given an eclipse, figure out the matching MCAT filename and then pull it from the bucket
    visit_metadata = pq.read_table(visit_metadata_fn)
    mcat_manifest = pd.read_csv(mcat_manifest_fn)
    row = visit_metadata.filter(pa.compute.equal(visit_metadata['ECLIPSE'], eclipse)).to_pydict()
    img, tilenum = int(row['IMG'][0]), int(row['TILENUM'][0])
    fp = mcat_manifest[mcat_manifest['source']=='visitI'][mcat_manifest['tilenum']==tilenum][mcat_manifest['img']==img]['fileNPath'].values[0]
    mcat_fn = fp.split('/')[-1]
    datapath = f"{rootpath}/e{str(eclipse).zfill(5)}"
    outpath = f"{datapath}/{mcat_fn}"
    if not os.path.exists(outpath):
        Bucket(bucket_name).get(mcat_fn,outpath)
    return outpath

def write_mcat_pos(eclipse,mcat_path,datapath):
    mcat_tbl = pdr.read(mcat_path)
    mcat_pos = mcat_tbl['GALEX_MERGED_SOURCE_LIST'][['alpha_j2000','delta_j2000']].rename(
        columns={'alpha_j2000':'ra','delta_j2000':'dec'})
    mcat_pos['eclipse']=np.full(len(mcat_pos),eclipse)
    mcat_pos_path = f'{datapath}/e{str(eclipse).zfill(5)}/e{str(eclipse).zfill(5)}_mcat_pos.csv'
    mcat_pos.to_csv(mcat_pos_path,index=None)
    return mcat_pos_path

def run_pipeline_from_cat(eclipse,mcat_cat_file,datapath,band='NUV',
                          aperture_sizes=[1.5, 2.3, 3.8, 6.0, 9.0, 12.8, 17.3],
                          compression="rice"):
    execute_pipeline(
        eclipse,
        band,
        depth=120,
        threads=4,
        local_root=datapath,
        recreate=False,
        aperture_sizes=aperture_sizes,
        write={"movie": True, "image": True},
        #coregister_lightcurves=True,
        photometry_only=False,
        compression=compression,
        suffix='mon',
        source_catalog_file=mcat_cat_file,
    )

def counts2mag(cps, band):
    scale = 18.82 if band == 'FUV' else 20.08
    with np.errstate(invalid='ignore'):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mag = -2.5 * np.log10(cps) + scale
    return mag

def mag2counts(mag, band):
    scale = 18.82 if band == 'FUV' else 20.08
    return 10.**(-(mag-scale)/2.5)

zpmag={'NUV':20.08, 'FUV':18.82}


