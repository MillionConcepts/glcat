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

    photonlist = savepath + photonlist_path

    if os.path.exists(photonlist):
        try:
            # get photonlist from bucket
            nf = parquet.read_table(photonlist,
                                    columns=['col', 'row', 'ra', 'dec', 't']).to_pandas()

            print(len(nf))

            nf['row_rnd'] = nf['row'].round().astype(int)
            nf['col_rnd'] = nf['col'].round().astype(int)

            # filtering photonlist for on detector
            nf = nf[(nf['col_rnd'] < 800) & (nf['row_rnd'] < 800)
                    & (nf['ra'] != 0) & (nf['dec'] != 0)
                    & (nf['col_rnd'] >= -50) & (nf['row_rnd'] >= -50)]
            mask = pd.notna(nf['ra'])
            nf = nf[mask]

            # rough approx not accounting for dead time
            expt = nf.iloc[len(nf) - 1]['t'] - nf.iloc[0]['t']

            # std dev of ra and dec by col and row w
            ra_stdev = bin2d(nf['col'], nf['row'], nf['ra'], 'std', nbins)
            dec_stdev = bin2d(nf['col'], nf['row'], nf['dec'], 'std', nbins)

            # photon count per bin
            count = bin2d(nf['col'], nf['row'], nf['ra'], 'count', nbins)
            count = count / expt

            # density mask
            density_mask = count >= .9
            # dispersion mask
            disp_mask = ra_stdev + dec_stdev > .014
            # response map
            dark_mask = count <= .015

            # empty masks
            hmask = np.ones(count.shape, dtype=bool)
            cmask = np.ones(count.shape, dtype=bool)

            # saving hotspot mask per eclipse and band
            hmask[density_mask & disp_mask] = 0
            hmask.tofile(f'{savepath}{eclipse}-{band}d-hmask.bin')

            # saving coldspot mask per eclipse and band
            cmask[dark_mask] = 0
            cmask.tofile(f'{savepath}{eclipse}-{band}d-cmask.bin')

        except KeyboardInterrupt:
            raise
        except Exception as ex:
            print(ex)
            print(f"failed {eclipse}")
    else:
        print("fail!")

    gc.collect()

    return
