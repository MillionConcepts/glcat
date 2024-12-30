if __name__ == "__main__":
    from backplanes import make_hotspot_files_eclipse_shorter
    import pandas as pd
    import gc
    nuv_list = pd.read_csv("every3_fuv.csv")

    # Access the column as a pandas Series object
    eclipses = nuv_list['eclipse'][300:]

    # Now, you can iterate through the values in the column_series
    for e in eclipses:
        e = int(e)
        make_hotspot_files_eclipse_shorter(
            eclipse=e,
            band="FUV"
        )
        gc.collect()