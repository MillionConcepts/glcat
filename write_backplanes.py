if __name__ == "__main__":
    from backplanes import make_backplanes

    make_backplanes(
        eclipse=21422,
        band="NUV",
        depth=1,
        leg=1,
        threads=4,
        burst=True,
        # should be root path for eclipse directories where you've written
        # your photonlists; will also write files
        # there. recall photonlists must be written with extended metadata
        local="/home/michael/Desktop/gphoton2/test_data",
        # "xy" or "dose"
        kind="xy",
        # detector radius cutoff -- only matters for dosemaps. pop up to ~750
        # if you want the stims.
        radius=400,
        # write arrays / xylists? (xylists relevant only for dosemaps)
        write={'array': True, 'xylist': True},
        # write arrays inline? possible only in burst mode.
        inline=True,
        # threshold and star_size for DAOStarFinder. matters only for
        # dosemap xylist
        threshold=0.75,
        star_size=2
    )

