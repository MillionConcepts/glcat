from pyarrow import parquet
import pandas as pd
import sys
import os
import shutil
from backplanes import make_backplanes
sys.path.append('/home/ubuntu/gPhoton2')
from gPhoton.pipeline import execute_pipeline


metadata = parquet.read_table('/home/ubuntu/gPhoton2/gPhoton/aspect/metadata.parquet',
                              columns=['eclipse', 'legs']).to_pandas()
band = "NUV"

fails = []

for eclipse in metadata['eclipse'][:10]:
    # how many legs does this eclipse have
    legs = metadata[metadata['eclipse'] == eclipse]['legs'].item()

    # gphoton run to produce extended photonlist
    # modified aspect table must be called "aspect2" and be in the aspect folder of gPhoton2
    try:
         execute_pipeline(
            eclipse,
            band,
            depth=120,
            # integer; None to deactivate (default None)
            threads=4,
            # where to both write output data and look for input data
            local_root="/home/ubuntu/gPhoton2/test_data",
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
    except KeyboardInterrupt:
        #raise
        print(f"keyboard interrupt :( for {eclipse} ")
        with open("failed_gphoton_eclipses.csv", "a+") as stream:
            stream.write(f"{eclipse},keyboard\n")

    except Exception as ex:
        print(f"something didn't work :( for {eclipse} ")
        with open("failed_gphoton_eclipses.csv", "a+") as stream:
            stream.write(f"{eclipse},{str(ex).replace(',', '')}\n")

    # if photonlist is produced, try to make backplanes
    try:
        legs = [0] if legs == 0 else tuple(range(legs))
        for leg in legs:
            make_backplanes(
                eclipse=eclipse,
                band=band,
                depth=1,
                leg=leg,
                threads=4,
                burst=True,
                local="/home/ubuntu/gPhoton2/test_data",
                kind="dose",
                radius=600,
                write={'array': True, 'xylist': False},
                inline=True,
                threshold=.45,
                star_size=2,
                snippet=None
                )
    except KeyboardInterrupt:
         #raise
        print(f"keyboard interrupt :( for {eclipse} ")
        with open("failed_backplane_eclipses.csv", "a+") as stream:
            stream.write(f"{eclipse},keyboard\n")
    except Exception as ex:
        print(f"something didn't work :( for {eclipse} ")
        with open("failed_backplane_eclipses.csv", "a+") as stream:
            stream.write(f"{eclipse},{str(ex).replace(',', '')}\n")

    pad_eclipse = str(eclipse).zfill(5)
    b = "n" if band == "NUV" else "f"

    #delete raw6
    #os.remove(f'home/ubuntu/gPhoton2/test_data/e{pad_eclipse}/e{pad_eclipse}-{b}d-raw6.fits.gz')

    # move folder of eclipse stuff (backplanes, extended photonlist) to s3
    try:
        dest = shutil.move(f'home/ubuntu/gPhoton2/test_data/e{pad_eclipse}', '/mnt/s3/')
        print(f"moved folder of {pad_eclipse} to {dest}")
    except KeyboardInterrupt:
        #raise
        print(f"keyboard interrupt :( for {eclipse} ")
        with open("failed_transfer_eclipses.csv", "a+") as stream:
            stream.write(f"{eclipse},keyboard\n")
    except Exception as ex:
        print("something didn't work :( ")
        with open("failed_transfer_eclipses.csv", "a+") as stream:
            stream.write(f"{eclipse},{str(ex).replace(',', '')}\n")

print("done :)")





