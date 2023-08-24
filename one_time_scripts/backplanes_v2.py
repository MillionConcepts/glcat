from pyarrow import parquet
import sys
import shutil
import subprocess
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
        local_root = "/home/ubuntu/gPhoton2/test_data",
        cmd = f"python pipeline_cli.py {eclipse} {band}" \
               f" --threads=4 --local_root={local_root} --stop_after='photonpipe'" \
               f" --compression=rice --write={{'movie':False,'image':False}} --extended_photonlist=True --aspect='aspect2' "
        subprocess.run(cmd, shell=True, timeout=300)
    except KeyboardInterrupt:
        raise
    except subprocess.TimeoutExpired as e:
        print(f"something didn't work :( for {eclipse} ")
        with open("failed_gphoton_eclipses.csv", "a+") as stream:
            stream.write(f"{eclipse},{str(e).replace(',', '')}\n")
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
        raise
    except Exception as ex:
        print(f"something didn't work :( for {eclipse} ")
        with open("failed_backplane_eclipses.csv", "a+") as stream:
            stream.write(f"{eclipse},{str(ex).replace(',', '')}\n")

    pad_eclipse = str(eclipse).zfill(5)
    b = "n" if band == "NUV" else "f"

    # move folder of eclipse stuff (backplanes, extended photonlist) to s3
    try:
        dest = shutil.move(f'home/ubuntu/gPhoton2/test_data/e{pad_eclipse}', '/mnt/s3/')
        print(f"moved folder of {pad_eclipse} to {dest}")
    except KeyboardInterrupt:
        raise
    except Exception as ex:
        print(f"something didn't work :( for {eclipse} ")
        with open("failed_transfer_eclipses.csv", "a+") as stream:
            stream.write(f"{eclipse},{str(ex).replace(',', '')}\n")

print("done :)")