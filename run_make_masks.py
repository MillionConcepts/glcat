if __name__ == "__main__":
    from backplanes import make_masks_per_eclipse
    import pandas as pd
    import gc
    import time

    paths = pd.read_csv("/glcat/notebooks/masks_flats/photonfile_keys.csv")
    eclipses = nuv_list['parq_keys'][0:]

    # band
    band = "n"

    savepath = '/mnt/s3/'

    nbins = 800 # works a little better than 800

    for photonlist_path in paths:
        start_time = time.time()

        eclipse = photonlist_path[0:6]

        make_masks_per_eclipse(
            eclipse,
            band,
            photonlist_path,
            nbins,
            savepath
        )
        end_time = time.time()

        duration = end_time - start_time
        print(f"that eclipse took {duration:.4f} seconds")

        gc.collect()

