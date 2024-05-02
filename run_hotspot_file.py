if __name__ == "__main__":
    from backplanes import make_hotspot_files_eclipse
    import pandas as pd

    nuv_list = pd.read_csv("every3_fuv.csv")

    # Access the column as a pandas Series object
    eclipses = nuv_list['eclipse']

    # Now, you can iterate through the values in the column_series
    for e in eclipses:
        make_hotspot_files_eclipse(
            eclipse=e,
            band="FUV"
        )

