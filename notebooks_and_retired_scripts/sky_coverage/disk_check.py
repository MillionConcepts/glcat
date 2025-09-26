"""Sanity check on the results of disks.py.  We should have three disks
in the disks files for every good aspect point."""

import argparse
import gzip
import multiprocessing
import struct
import sys
import time

from datetime import timedelta
from pathlib import Path

import numpy as np

from pyarrow import parquet


START_TIME = None
def progress(msg: str) -> None:
    now = time.monotonic()

    global START_TIME
    if START_TIME is None:
        START_TIME = now

    elapsed = timedelta(seconds=now - START_TIME)
    sys.stderr.write(f"[{elapsed}]  {msg}\n")


def count_aspect_fixes(aspect):
    progress("reading aspect")
    all_aspect = parquet.read_table(
        aspect,
        columns=["eclipse", "time", "flags", "ra", "dec"]
    ).to_pandas().groupby("eclipse")

    progress("counting good fixes")
    eclipses = 0
    fixes = 0
    for eclipse, aspect in all_aspect:
        ok = aspect["flags"].to_numpy() & 1 == 0

        # we also distrust aspect records immediately before and after a
        # flagged aspect record; the first record in 'aspect' is assumed
        # to be preceded by a flagged record, and the last record is
        # assumed to be followed by a flagged record
        ok &= np.concat([[False], ok[:-1]])
        ok &= np.concat([ok[1:], [False]])

        # we also distrust aspect records whose time step is not 1
        times = aspect["time"].to_numpy()
        ok[1:] &= np.isclose(times[1:] - times[:-1], 1.0)

        eclipses += 1
        fixes += ok.sum()

        if eclipses % 10000 == 0:
            progress(f"{eclipses}/34389 eclipses = {fixes} fixes")

    progress(f"total {fixes} fixes")
    return fixes


def count_recorded_disks_1(fname):
    fixes = 0
    errors = False
    with gzip.open(fname, "rb") as fp:
        try:
            while True:
                n = fp.read(4)
                if len(n) == 0:
                    break
                if len(n) < 4:
                    raise ValueError("truncated record (length field)")
                n = struct.unpack("<I", n)[0]
                lwkb = (n >> 16) & 0x3FFF
                data = fp.read(lwkb)
                if len(data) < lwkb:
                    raise ValueError("truncated record (payload)")
                fixes += 1
            progress(f"{fname.name}: {fixes} disk sets")
        except Exception as e:
            progress(f"{fname.name}: error: {e} after {fixes} disk sets")
            errors = True

    return fixes, errors

def count_recorded_disks(disksdir):
    progress("counting recorded disk sets")
    fixes = 0
    errors = False
    with multiprocessing.Pool() as pool:
        for f_fixes, f_errors in pool.imap_unordered(
                count_recorded_disks_1,
                disksdir.glob("disks-*.gz"),
                chunksize=1,
        ):
            fixes += f_fixes
            errors |= f_errors
    progress(f"total {fixes} disk sets")
    return fixes, errors


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("aspect", type=Path)
    ap.add_argument("outdir", type=Path)
    args = ap.parse_args()
    aspect = args.aspect
    outdir = args.outdir

    expected_fixes = count_aspect_fixes(aspect)
    got_fixes, errors = count_recorded_disks(outdir)
    if expected_fixes == got_fixes:
        sys.stderr.write(f"ok! expected = got = {expected_fixes}\n")
    else:
        sys.stderr.write(
            f"*** error: expected {expected_fixes}, got {got_fixes}\n"
        )
        errors = True

    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()

