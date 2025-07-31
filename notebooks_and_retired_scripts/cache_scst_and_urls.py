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
from limiter import Limiter
from pyarrow import parquet, Table, Array, ChunkedArray


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
        scst_file = mast.download_data(eclipse, ftype='scst', datadir=scst_dir)
        return MastMetadata(
            eclipse = eclipse,
            scst_file = scst_file,
            nuv_raw6 = paths.get("NUV") is not None,
            fuv_raw6 = paths.get("FUV") is not None,
        )


def eclipses_to_do(
    datadir: Path,
    shelf: shelve.Shelf,
    metadata: Table
) -> set[int]:
    eclipses = set(e.as_py() for e in metadata.column('eclipse'))
    prune = []
    for e, mm in shelf.items():
        try:
            scst = datadir / mm["scst"]
            st = scst.stat()
            if st.st_size > 0:
                eclipses.remove(int(e))
                continue

            # if we reach this point, the scst file exists but is truncated;
            # delete it so mast.download_data will retry
            scst.unlink()
        except (KeyError, FileNotFoundError):
            pass

        # if we reach this point, there's something wrong with the
        # record in the shelf, and we should regenerate it
        prune.append(e)

    for e in prune:
        del shelf[e]

    return eclipses


def process_eclipses(
    datadir: Path,
    metadata: Table,
    average_rate: int,
    burst_rate: int,
) -> None:

    progress("counting eclipses...")
    rate_limit = Limiter(average_rate, burst_rate)
    datadir.mkdir(parents=True, exist_ok=True)

    with ThreadPoolExecutor() as pool, \
         shelve.open(datadir / "index.shelf") as shelf:
        eclipses = eclipses_to_do(datadir, shelf, metadata)

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
                shelf[str(mm.eclipse)] = {
                    "scst": str(mm.scst_file.relative_to(datadir)),
                    "nuv_raw6": mm.nuv_raw6,
                    "fuv_raw6": mm.fuv_raw6
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
    ap.add_argument("metadata", type=Path,
                    help="Parquet file containing eclipse metadata")
    ap.add_argument("data_dir", type=Path,
                    help="Where to write SCST files and index")
    ap.add_argument("--avg-rate", type=int, default=500,
                    help="Average number of MAST queries per second")
    ap.add_argument("--burst-rate", type=int, default=1000,
                    help="Burst rate (absolutely no more than this many queries per second)")
    args = ap.parse_args()

    progress("loading metadata table...")
    md = parquet.read_table(args.metadata)
    if "eclipse" not in md.column_names:
        sys.stderr.write(f"{args.metadata}: error: no 'eclipse' column found")
        sys.exit(1)

    process_eclipses(args.data_dir, md, args.avg_rate, args.burst_rate)


main()
