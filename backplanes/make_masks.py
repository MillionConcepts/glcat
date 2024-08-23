import sys
import os
import shutil
import gc
from quickbin import bin2d
import pandas as pd
import numpy as np
import pyarrow.parquet as parquet


def make_masks_per_eclipse(eclipse, band, photonlist_path, nbins, savepath):
    # use existing photonlists to make individual hotspot and coldspot masks
    # per eclipse band

    photonlist = str(savepath + photonlist_path)
    print(photonlist)
    if os.path.exists(photonlist):
        try:
            # get photonlist from bucket
            print("reading photonlist")

            parquet_file = parquet.ParquetFile(photonlist)
            nrows = parquet_file.metadata.num_rows
            if nrows > 20000000:
                n = int(nrows / 20000000)
            else:
                n = 1
            nf = pd.DataFrame()

            for chunk in filter_parquet_with_iter_batches(photonlist, n):
                print("adding chunk")
                nf = pd.concat([nf, chunk])

            #nf = parquet.read_table(photonlist,
            #                        columns=['col', 'row', 'ra', 'dec', 't']).to_pandas()\
            #for chunk in pd.read_parquet(photonlist,
            #                             columns=['col', 'row', 'ra', 'dec', 't'],
            #                             chunksize=n):
           #     nf = pd.concat([nf, chunk.iloc[::0]])

            nf['row_rnd'] = nf['row'].round().astype(int)
            nf['col_rnd'] = nf['col'].round().astype(int)

            print("filtering photonlist")
            # filtering photonlist for on detector
            nf = nf[(nf['col_rnd'] < 800) & (nf['row_rnd'] < 800)
                    & (nf['ra'] != 0) & (nf['dec'] != 0)
                    & (nf['col_rnd'] >= -50) & (nf['row_rnd'] >= -50)]
            mask = pd.notna(nf['ra'])
            nf = nf[mask]

            print("calculating expt")
            # rough approx not accounting for dead time
            expt = nf.iloc[len(nf) - 1]['t'] - nf.iloc[0]['t']/n

            print("quickbinning")
            # std dev of ra and dec by col and row
            ra_stdev = bin2d(nf['col'], nf['row'], nf['ra'], 'std', nbins)
            dec_stdev = bin2d(nf['col'], nf['row'], nf['dec'], 'std', nbins)

            # photon count per bin
            count = bin2d(nf['col'], nf['row'], nf['ra'], 'count', nbins)
            count = count / expt

            print("masking binned data")
            # density mask
            density_mask = count >= .9
            # dispersion mask
            disp_mask = ra_stdev + dec_stdev > .014
            # response map
            dark_mask = count <= .015

            print("making new masks")
            # empty masks
            hmask = np.ones(count.shape, dtype=bool)
            cmask = np.ones(count.shape, dtype=bool)

            # saving hotspot mask per eclipse and band
            hmask[density_mask & disp_mask] = 0
            print("saving hotspot mask")
            hmask.tofile(f'{savepath}{eclipse}-{band}d-hmask.bin')

            # saving coldspot mask per eclipse and band
            cmask[dark_mask] = 0
            print("saving coldspot mask")
            cmask.tofile(f'{savepath}{eclipse}-{band}d-cmask.bin')

        except KeyboardInterrupt:
            raise
        except Exception as ex:
            print(ex)
            print(f"failed {eclipse}")
    else:
        print("fail!")

    print("cleaning up")
    gc.collect()

    return


def filter_parquet_with_iter_batches(file_path, batch_size):
    parquet_file = parquet.ParquetFile(file_path)
    for batch in parquet_file.iter_batches(columns=['col', 'row', 'ra', 'dec', 't'],
                                           batch_size=batch_size,
                                           filters=[('col', '<=', 800),
                                                    ('row', '<=', 800),
                                                    ('col', '>=', -50),
                                                    ('row', '>=', -50)]):
        df = batch.to_pandas()
        filtered_df = df.iloc[[0]] if not df.empty else pd.DataFrame()
        yield filtered_df
