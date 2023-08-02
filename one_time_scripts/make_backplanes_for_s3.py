from pyarrow import parquet
import pandas as pd
import subprocess
from glcat.backplanes import make_backplanes
from gPhoton2.pipeline import execute_pipeline

metadata = parquet.read_table('/gPhoton2/aspect/metadata.parquet',
                              columns=['eclipse', 'legs']).to_pandas()
band = "NUV"

fails = []

for eclipse in metadata['eclipse']:
    # how many legs does this eclipse have
    legs = metadata[metadata['eclipse'] == eclipse]['legs'].item()

    # gphoton run to produce extended photonlist
    # modified aspect table must be called "aspect2" and be in the aspect folder of gPhoton2
    execute_pipeline(
            eclipse,
            band,
            depth=120,
            # integer; None to deactivate (default None)
            threads=4,
            # where to both write output data and look for input data
            local_root="test_data",
            # auxiliary remote location for input data
            # remote_root="/mnt/s3",
            recreate=True,
            # list of floats; relevant only to lightcurve / photometry portion
            aperture_sizes=[12.8],
            # actually write image/movie products? otherwise hold in memory but
            # discard (possibly after performing photometry).
            write={"movie": False, "image": False},
            coregister_lightcurves=False,
            # photonpipe, moviemaker, None (default None)
            stop_after='photonpipe',
            photometry_only=False,
            # None, "gzip", "rice"
            compression="rice",
            # use array sparsification on movie frames?
            lil=True,
            # write movie frames as separate files
            burst=False,
            extended_photonlist=True,
            # aspect file, don't need to set unless need to use alt
            # file, 'aspect2.parquet'
            aspect="aspect2"
        )

    # if photonlist is produced, try to make backplanes
    try:
        if legs > 0:
            for leg in range(legs):
                make_backplanes(
                        eclipse=eclipse,
                        band=band,
                        depth=1,
                        leg=leg,
                        threads=4,
                        burst=True,
                        local="/gPhoton2",
                        kind="dose",
                        radius=600,
                        write={'array': True, 'xylist': False},
                        inline=True,
                        threshold=.45,
                        star_size=2,
                        snippet=None
                )
        else:
            make_backplanes(
                eclipse=eclipse,
                band=band,
                depth=1,
                leg=0,
                threads=4,
                burst=True,
                local="/gPhoton2",
                kind="dose",
                radius=600,
                write={'array': True, 'xylist': False},
                inline=True,
                threshold=.45,
                star_size=2,
                snippet=None
                )
    except:
        print("something didn't work :( ")
        fails = fails.append(eclipse)

    pad_eclipse = str(eclipse).zfill(5)
    # move folder of eclipse stuff (backplanes, extended photonlist) to s3
    move = f"mv e{pad_eclipse} /mnt/s3/backplanes/"
    subprocess.call(move, shell=True)

# save failures
fails_df = pd.DataFrame(fails)
fails_df.to_csv("failed_backplane_eclipses.csv")





