import sys
import os
import shutil
from backplanes import make_backplanes
sys.path.append('/home/ubuntu/gPhoton2')
from gPhoton.pipeline import execute_pipeline
import gc


def make_hotspot_files_eclipse_shorter(eclipse, band):
    # gphoton run to produce extended photonlist
    # modified aspect table must be called "aspect2" and be in the aspect folder of gPhoton2

    # try to run gphoton
    pad_eclipse = str(eclipse).zfill(5)
    b = "n" if band == "NUV" else "f"

    photonlist_path = f"/mnt/s3/e{pad_eclipse}-{b}d--NF/e{pad_eclipse}-{b}d-b00.parquet"

    if os.path.exists(photonlist_path):
        print(f"File {photonlist_path} exists.")
        # then try to make full depth backplane
        try:
            make_backplanes(
                eclipse=eclipse,
                band=band,
                depth=100,
                leg=0,
                threads=4,
                burst=True,
                local=f'/mnt/s3/e{pad_eclipse}-{b}d--NF',
                kind="dose",
                radius=600,
                write={'array': True, 'xylist': False},
                inline=True,
                threshold=.45,
                star_size=2,
                snippet=None)
            gc.collect()

            #try:
            #    dest = shutil.move(f'/home/ubuntu/gPhoton2/test_data/e{pad_eclipse}', f'/mnt/s3/e{pad_eclipse}-{b}d--NF')
            #    print(f"moved folder of {pad_eclipse} to {dest}")
        except KeyboardInterrupt:
            raise
        except Exception as ex:
            print(f"failed {eclipse}, {ex}")

    if os.path.exists(f'/home/ubuntu/gPhoton2/test_data/temp'):
        print("deleting temp folder from gphoton test data")
        shutil.rmtree(f'/home/ubuntu/gPhoton2/test_data/temp')

    if os.path.exists(f'/home/ubuntu/gPhoton2/test_data/e{pad_eclipse}'):
        print("deleting failed folder from ec2")
        shutil.rmtree(f'/home/ubuntu/gPhoton2/test_data/e{pad_eclipse}')

    return
