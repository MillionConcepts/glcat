import argparse
import sys
import time

from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import timedelta
from enum import Flag, auto
from typing import Iterable

import numpy as np

from gPhoton.io.mast import get_raw_paths
from limiter import Limiter
from pyarrow import parquet, Table


START_TIME = None
def progress(msg):
    now = time.monotonic()

    global START_TIME
    if START_TIME is None:
        START_TIME = now

    elapsed = timedelta(seconds=now - START_TIME)
    sys.stderr.write(f"[{elapsed}]  {msg}\n")


class Band(Flag):
    NONE = 0
    NUV  = auto()
    FUV  = auto()
    BOTH = NUV | FUV


def has_band(eband: Band | None, band: Band) -> bool | None:
    """True if 'band' is included in 'eband', False if it isn't,
       and None if eband is None."""
    if eband is None: return None
    return (eband & band) == band

def available_bands(eclipse: int, rate_limit: Limiter) -> (int, Band):
    """Ask MAST whether data was collected in each of the detector
       bands for the given eclipse."""

    with rate_limit:
        paths = get_raw_paths(eclipse)

    bands = Band.NONE
    if paths.get("NUV") is not None:
        bands |= Band.NUV
    if paths.get("FUV") is not None:
        bands |= Band.FUV
    return (eclipse, bands)


def bands_for_eclipses(
        eclipses: Iterable[int],
        *,
        average_rate: int,
        burst_rate: int,
) -> dict[int, Band]:

    rate_limit = Limiter(average_rate, burst_rate)
    with ThreadPoolExecutor() as pool:

        futures = []
        progress("distributing eclipses...")
        for eclipse in eclipses:
            if len(futures) > 0 and len(futures) % 1000 == 0:
                progress(f"{len(futures)} eclipses distributed")
            futures.append(pool.submit(available_bands, eclipse, rate_limit))
        if len(futures) > 0 and len(futures) % 1000 != 0:
            progress(f"{len(futures)} eclipses distributed")

        progress("querying eclipses...")
        bands = {}
        errors = 0
        for fut in as_completed(futures):
            if len(bands) > 0 and len(bands) % 1000 == 0:
                progress(f"{len(bands)} eclipses queried, {errors} errors")
            try:
                e, b = fut.result()
                bands[e] = b
            except Exception as e:
                t = type(e)
                m = str(e)
                if t == m:
                    sys.stderr.write(f"{eclipse}: error: {m}\n")
                else:
                    sys.stderr.write(f"{eclipse}: error: {t}: {m}\n")
                bands[eclipse] = None
                errors += 1

        if len(bands) > 0 and len(bands) % 1000 != 0:
            progress(f"{len(bands)} eclipses queried, {errors} errors")
        return bands


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("metadata", help="Parquet file containing eclipse metadata")
    ap.add_argument("output", help="Where to write augmented metadata file")
    ap.add_argument("--avg-rate", type=int, default=30,
                    help="Average number of MAST queries per second")
    ap.add_argument("--burst-rate", type=int, default=100,
                    help="Burst rate (absolutely no more than this many queries per second)")
    args = ap.parse_args()

    progress("loading metadata table")
    md_in = parquet.read_table(args.metadata)
    if "eclipse" not in md_in.column_names:
        sys.stderr.write(f"{args.metadata}: error: no 'eclipse' column found")
        sys.exit(1)

    progress("preprocessing metadata table")
    if "have_nuv" in md_in.column_names:
        md_in = md_in.drop_columns(["have_nuv"])
    if "have_fuv" in md_in.column_names:
        md_in = md_in.drop_columns(["have_fuv"])

    md_in = md_in.sort_by("eclipse")
    eclipses = md_in.column("eclipse").to_numpy()

    progress(f"{len(eclipses)} eclipses to process")
    eclipse_bands = bands_for_eclipses(
        eclipses, average_rate=args.avg_rate, burst_rate=args.burst_rate
    )

    progress("making column vectors")
    have_nuv = [
        has_band(eclipse_bands.get(eclipse), Band.NUV)
        for eclipse in eclipses
    ]
    have_fuv = [
        has_band(eclipse_bands.get(eclipse), Band.FUV)
        for eclipse in eclipses
    ]
    update = Table.from_pydict({
        "eclipse": eclipses,
        "have_nuv": have_nuv,
        "have_fuv": have_fuv
    })

    progress("merging tables")
    md_out = md_in.join(update, keys=["eclipse"], join_type="left outer")
    md_out = md_out.sort_by("eclipse")

    progress("writing output")
    parquet.write_table(md_out, args.output)


main()
