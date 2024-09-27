if __name__ == "__main__":
    from backplanes import make_masks_per_eclipse
    import pandas as pd
    import gc
    import time
    import pyarrow.parquet as parquet

    paths = pd.read_csv("~/glcat/notebooks/masks_flats/photonfile_keys.csv")
    eclipses = paths['parq_keys'][0:]

    cutoff_data = parquet.read_table('~/glcat/notebooks/masks_flats/stdev_by_eclipse.parquet').to_pandas()

    # band
    band = "n"

    savepath = '/mnt/s3/'

    nbins = 800 # works a little better than 800

    for photonlist_path in eclipses:
        print(f"running {photonlist_path}")
        start_time = time.time()

        eclipse = photonlist_path[0:6]

        ra_cutoff = cutoff_data[(cutoff_data['eclipse'] == eclipse)]['ra']
        dec_cutoff = cutoff_data[(cutoff_data['eclipse'] == eclipse)]['dec']
        print(f"cutoffs: {ra_cutoff}, {dec_cutoff}")
        make_masks_per_eclipse(
            eclipse,
            band,
            ra_cutoff,
            dec_cutoff,
            photonlist_path,
            nbins,
            savepath
        )
        end_time = time.time()

        duration = end_time - start_time
        print(f"that eclipse took {duration:.4f} seconds")

        gc.collect()

