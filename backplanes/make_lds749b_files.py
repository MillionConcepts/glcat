import sys
import os
import shutil
from backplanes import make_backplanes
sys.path.append('/home/ubuntu/gPhoton2')
from gPhoton.pipeline import execute_pipeline
import gc


def gphoton_only(eclipse, band):
    # gphoton run to produce extended photonlist

    # try to run gphoton
    pad_eclipse = str(eclipse).zfill(5)
    b = "n" if band == "NUV" else "f"


    try:
        print(f"running gphoton2 on eclipse {pad_eclipse} with new mask")
        execute_pipeline(
            eclipse,
            band,
            depth=None,
            # integer; None to deactivate (default None)
            threads=None,
            # where to both write output data and look for input data
            local_root="/home/ubuntu/gPhoton2/test_data",
            # auxiliary remote location for input data
            # remote_root="/mnt/s3",
            recreate=True,
            # list of floats; relevant only to lightcurve / photometry portion
            aperture_sizes=[12.8],
            # actually write image/movie products? otherwise hold in memory but
            # discard (possibly after performing photometry).
            write={"movie": False, "image": True},
            coregister_lightcurves=False,
            # photonpipe, moviemaker, None (default None)
            stop_after=None,
            photometry_only=False,
            # None, "gzip", "rice"
            compression="rice",
            # use array sparsification on movie frames?
            lil=True,
            # write movie frames as separate files
            burst=False,
            extended_photonlist=True)
        gc.collect()

        # then try to make full depth backplane
        try:
            make_backplanes(
                eclipse=eclipse,
                band=band,
                depth=1000,
                leg=0,
                threads=4,
                burst=True,
                local="/home/ubuntu/gPhoton2/test_data",
                kind="dose",
                radius=600,
                write={'array': True, 'xylist': False},
                inline=True,
                threshold=.45,
                star_size=2,
                snippet=None)
            gc.collect()

        except KeyboardInterrupt:
            raise
        except Exception as ex:
            print(f"failed backplane {eclipse} ")
            with open("failed_backplane_eclipses.csv", "a+") as stream:
                stream.write(f"{eclipse},{str(ex).replace(',', '')}\n")

        try:
            dest = shutil.move(f'/home/ubuntu/gPhoton2/test_data/e{pad_eclipse}', f'/mnt/s3/e{pad_eclipse}-{b}d-lds749b')
            print(f"moved folder of {pad_eclipse} to {dest}")
        except KeyboardInterrupt:
            raise
        except Exception as ex:
            print(f"failed transfer {eclipse} ")
            with open("failed_transfer_eclipses.csv", "a+") as stream:
                stream.write(f"{eclipse},{str(ex).replace(',', '')}\n")

    except KeyboardInterrupt:
        raise
    except Exception as ex:
        print(f"failed gphoton {eclipse} ")
        with open("failed_gphoton_eclipses.csv", "a+") as stream:
            stream.write(f"{eclipse},{str(ex).replace(',', '')}\n")

    if os.path.exists(f'/home/ubuntu/gPhoton2/test_data/temp'):
        print("deleting temp folder from gphoton test data")
        shutil.rmtree(f'/home/ubuntu/gPhoton2/test_data/temp')

    if os.path.exists(f'/home/ubuntu/gPhoton2/test_data/e{pad_eclipse}'):
        print("deleting failed folder from ec2")
        shutil.rmtree(f'/home/ubuntu/gPhoton2/test_data/e{pad_eclipse}')

    return
