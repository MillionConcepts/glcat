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
    # per eclipse & band

    photon_file = str(savepath + photonlist_path)
    print(photon_file)
    if os.path.exists(photon_file):
        try:
            # get photonlist from bucket
            print("reading photonlist")
            photonlist = parquet.ParquetFile(photon_file)
            if photonlist.metadata.num_rows < 40000000: # kind of an arbitrary cutoff

                nf = photonlist.read_row_groups([0],columns=['col', 'row', 'ra', 'dec', 't']).to_pandas()
                print(len(nf))

                nf['row_rnd'] = nf['row'].round().astype(int)
                nf['col_rnd'] = nf['col'].round().astype(int)

                print("filtering photonlist")
                # filtering photonlist for on detector because read_row_groups can't
                nf = nf[(nf['col_rnd'] <= 800) & (nf['row_rnd'] <= 800)
                        & (nf['ra'] != 0) & (nf['dec'] != 0)
                        & (nf['col_rnd'] >= 0) & (nf['row_rnd'] >= 0)]
                mask = pd.notna(nf['ra'])
                nf = nf[mask]

                print("calculating expt")
                # rough approx not accounting for dead time
                expt = nf.iloc[len(nf) - 1]['t'] - nf.iloc[0]['t']

                ra = nf.iloc[0]['ra']
                dec = nf.iloc[0]['dec']
                t = nf.iloc[0]['t']

                # adding edge points
                edge_points = pd.DataFrame({
                    'col': [0, 0, 800, 800],
                    'row': [0, 800, 0, 800],
                    'ra': [ra,ra,ra,ra],
                    'dec': [dec,dec,dec,dec],
                    't': [t,t,t,t]
                })

                # Append the new rows to the DataFrame
                nf = nf.append(edge_points, ignore_index=True)

                print("quickbinning")
                ra_stdev = bin2d(nf['col'], nf['row'], nf['ra'], 'std', nbins)
                dec_stdev = bin2d(nf['col'], nf['row'], nf['dec'], 'std', nbins)
                count = bin2d(nf['col'], nf['row'], nf['ra'], 'count', nbins)
                count = count / expt

                print("masking binned data")
                density_mask = count >= .9
                disp_mask = ra_stdev + dec_stdev > .014
                dark_mask = count <= .008

                print("making new masks")
                hmask = np.ones(count.shape, dtype=bool)
                cmask = np.ones(count.shape, dtype=bool)

                hmask[density_mask & disp_mask] = 0
                print("saving hotspot mask")
                hmask.tofile(f'{savepath}{eclipse}-{band}d-hmask.bin')

                cmask[dark_mask] = 0
                print("saving coldspot mask")
                cmask.tofile(f'{savepath}{eclipse}-{band}d-cmask.bin')
            else:
                print("photonlist too big")
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

