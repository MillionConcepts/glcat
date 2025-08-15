import argparse
import sys
import shelve
import time

from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import timedelta
from pathlib import Path
from typing import Iterable

from gPhoton.io import mast
from astropy.io import fits
from limiter import Limiter
from pyarrow import parquet, Table, Array, ChunkedArray

DEFAULT_MAX_ECLIPSE = 50000

START_TIME = None
def progress(msg: str) -> None:
    now = time.monotonic()

    global START_TIME
    if START_TIME is None:
        START_TIME = now

    elapsed = timedelta(seconds=now - START_TIME)
    sys.stderr.write(f"[{elapsed}]  {msg}\n")


@dataclass
class MastMetadata:
    eclipse: int
    scst_file: Path
    nuv_raw6: bool
    fuv_raw6: bool


def retrieve_urls_and_scst(
    eclipse: int,
    scst_dir: Path,
    rate_limit: Limiter
) -> MastMetadata:
    with rate_limit:
        paths = mast.get_raw_paths(eclipse)
        if paths.get("scst") is None:
            scst_file = None
        else:
            scst_file = mast.download_data(
                eclipse, ftype='scst', datadir=scst_dir
            )
            try:
                with fits.open(scst_file, lazy_load_hdus=False) as hdulist:
                    pass
            except Exception:
                # the scst file is corrupt on MAST, rename it out of
                # the way and mark this eclipse as not having an scst file
                scst_file.rename(
                    scst_file.with_name(
                        scst_file.name + ".corrupt"
                    )
                )
                scst_file = None

        return MastMetadata(
            eclipse = eclipse,
            scst_file = scst_file,
            nuv_raw6 = paths.get("NUV") is not None,
            fuv_raw6 = paths.get("FUV") is not None,
        )


def eclipses_to_do(
    datadir: Path,
    shelf: shelve.Shelf,
    max_eclipse: int,
    eclipses_with_aspect: frozenset,
) -> set[int]:

    todo = set()
    for e in range(max_eclipse):
        se = str(e)
        if se not in shelf:
            todo.add(e)
            continue

        mm = shelf[se]
        if (
                "scst" not in mm
                or "nuv_raw6" not in mm
                or "fuv_raw6" not in mm
        ):
            # there's something wrong with the record in the shelf
            scst_1 = f"e{e:05}-scst.fits.gz"
            scst_2 = mm.get("scst")
            (datadir / scst_1).unlink()
            if scst_2 != scst_1 and scst_2 is not None:
                (datadir / scst_2).unlink()
            del shelf[se]
            todo.add(e)
            continue

        if mm["scst"] is not None:
            scst = datadir / mm["scst"]
            try:
                st = scst.stat()
            except FileNotFoundError:
                del shelf[se]
                todo.add(e)
                continue
            if st.st_size == 0:
                # the scst file exists but is truncated;
                # delete it so mast.download_data will retry
                scst.unlink()
                del shelf[se]
                todo.add(e)
                continue

            try:
                with fits.open(scst, lazy_load_hdus=False) as hdulist:
                    pass
            except Exception:
                # the scst file is corrupt, delete it so mast.download_data
                # will retry
                scst.unlink()
                del shelf[se]
                todo.add(e)
                continue

        if "gphoton_aspect" not in mm:
            mm["gphoton_aspect"] = (e in eclipses_with_aspect)
            shelf[se] = mm

    return todo


def process_eclipses(
    datadir: Path,
    average_rate: int,
    burst_rate: int,
    max_eclipse: int,
    eclipses_with_aspect: frozenset
) -> None:

    progress("counting eclipses...")
    rate_limit = Limiter(average_rate, burst_rate)
    datadir.mkdir(parents=True, exist_ok=True)

    with ThreadPoolExecutor() as pool, \
         shelve.open(datadir / "index.shelf") as shelf:
        eclipses = eclipses_to_do(
            datadir,
            shelf,
            max_eclipse,
            eclipses_with_aspect,
        )

        progress("queuing...")
        futures = []
        progress(f"0/{len(eclipses)} eclipses queued")
        for eclipse in eclipses:
            if len(futures) > 0 and len(futures) % 1000 == 0:
                progress(f"{len(futures)}/{len(eclipses)} eclipses queued")
            futures.append(pool.submit(
                retrieve_urls_and_scst,
                eclipse, datadir, rate_limit
            ))
        if len(futures) > 0 and len(futures) % 1000 != 0:
            progress(f"{len(futures)}/{len(eclipses)} eclipses queued")

        progress("querying...")
        completed = 0
        errors = 0
        for fut in as_completed(futures):
            if completed > 0 and completed % 1000 == 0:
                progress(f"{completed}/{len(futures)} eclipses queried, {errors} errors")
            try:
                mm = fut.result()
                if mm.scst_file is None:
                    scst = None
                else:
                    scst = str(mm.scst_file.relative_to(datadir))
                shelf[str(mm.eclipse)] = {
                    "scst": scst,
                    "nuv_raw6": mm.nuv_raw6,
                    "fuv_raw6": mm.fuv_raw6,
                    "gphoton_aspect": (mm.eclipse in eclipses_with_aspect),
                }
            except Exception as e:
                t = type(e)
                m = str(e)
                if t == m:
                    sys.stderr.write(f"{eclipse}: error: {m}\n")
                else:
                    sys.stderr.write(f"{eclipse}: error: {t}: {m}\n")
                errors += 1
            completed += 1

        if completed > 0 and completed % 1000 != 0:
            progress(f"{completed}/{len(futures)} eclipses queried, {errors} errors")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("aspect", type=Path,
                    help="Parquet file containing eclipse aspect data")
    ap.add_argument("data_dir", type=Path,
                    help="Where to write SCST files and index")
    ap.add_argument("--max-eclipse", type=int, default=DEFAULT_MAX_ECLIPSE,
                    help="Maximum eclipse number to query")
    ap.add_argument("--avg-rate", type=int, default=500,
                    help="Average number of MAST queries per second")
    ap.add_argument("--burst-rate", type=int, default=1000,
                    help="Burst rate (absolutely no more than this many queries per second)")
    args = ap.parse_args()

    progress("loading aspect table...")
    aspect = parquet.read_table(args.aspect)
    if "eclipse" not in aspect.column_names:
        sys.stderr.write(f"{args.metadata}: error: no 'eclipse' column found")
        sys.exit(1)

    eclipses_with_aspect = frozenset(
        aspect.column("eclipse").to_pylist()
    )

    process_eclipses(
        args.data_dir,
        args.avg_rate,
        args.burst_rate,
        args.max_eclipse,
        eclipses_with_aspect,
    )


main()
