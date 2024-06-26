if __name__ == "__main__":
    from backplanes import gphoton_only
    import pandas as pd
    import gc


    elist = pd.read_csv("lds749b.csv")

    # Access the column as a pandas Series object
    eclipses = elist['eclipse'][0:]

    # Now, you can iterate through the values in the column_series
    for e in eclipses:
        e = int(e)
        gphoton_only(
            eclipse=e,
            band="FUV"
        )
        gc.collect()
