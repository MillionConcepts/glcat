if __name__ == "__main__":
    from backplanes import make_backplanes

    make_backplanes(
        eclipse=23456,
        band="NUV",
        depth=120,
        leg=0,
        threads=4,
        burst=False,
        # should be root path for eclipse directories where you've written
        # your photonlists; will also write files
        # there. recall photonlists must be written with extended metadata
        local="/home/michael/Desktop/gphoton2/test_data",
        # "xy" or "dose"
        kind="dose",
        # detector radius cutoff -- only matters for dosemaps. pop up to ~750
        # if you want the stims.
        radius=400
    )

