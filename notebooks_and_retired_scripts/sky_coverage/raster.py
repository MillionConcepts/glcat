"""Phase 3 of sky coverage map generation.

For each octant of the celestial sphere, read in the point set for that
octant, then count the number of disks that overlap each point.  Maintain
six separate counts: three sizes of disk * two detector bands.
"""

import argparse
import collections
import gzip
import multiprocessing
import os
import pickle
import struct
import sys
import time

from contextlib import closing
from datetime import timedelta
from pathlib import Path

import numpy as np
import pyarrow as pa
import shapely

from pyarrow import parquet

RASTER_SCHEMA = pa.schema([
    ("gal_l", pa.float64()),
    ("gal_b", pa.float64()),
    ("nuv_350", pa.uint32()),
    ("nuv_370", pa.uint32()),
    ("nuv_400", pa.uint32()),
    ("fuv_350", pa.uint32()),
    ("fuv_370", pa.uint32()),
    ("fuv_400", pa.uint32()),
])


START_TIME = None
def progress(msg: str) -> None:
    now = time.monotonic()

    global START_TIME
    if START_TIME is None:
        START_TIME = now

    elapsed = timedelta(seconds=now - START_TIME)
    sys.stderr.write(f"[{elapsed}]  {msg}\n")


def disks_stream(fname):
    """Yield (eclipse, has_nuv, has_fuv, disks) for each disk set in
       the file FNAME.  disks = [disk_350, disk_370, disk_400].
    """
    with gzip.open(fname, "rb") as fp:
        try:
            while True:
                eclipse, has_nuv, has_fuv, wkbs = pickle.load(fp)
                yield (
                    eclipse, has_nuv, has_fuv,
                    [shapely.from_wkb(wkb) for wkb in wkbs]
                )
        except EOFError:
            pass


def flush_hits(hits, counts, drain):
    """
    Flush accumulated hits from HITS into COUNTS, avoiding taking locks
    unless either there is enough work to justify it, or we are being
    asked to drain the entire hits record.
    """
    for band, bhits in hits.items():
        for dsize, dhits in enumerate(bhits):
            if len(dhits) >= 5000 or drain:
                dcounts = counts[band][dsize]
                with dcounts.get_lock():
                    v = dcounts.get_obj()
                    for h, n in dhits.items():
                        v[h] += n
                dhits.clear()


def crunch_disks(worker_id, fname, pointset, counts):
    """
    Load disks from FNAME.  For each disk, find the members of the point
    set that are within it.  Count the times that each member of the point
    set was within some disk.
    """
    ctr = collections.Counter
    hits = {
        "n": [ ctr(), ctr(), ctr() ],
        "f": [ ctr(), ctr(), ctr() ],
    }
    prev_eclipse = None
    eclipses = 0
    hitcount = 0
    last_flush_hitcount = 0
    for eclipse, has_nuv, has_fuv, disks in disks_stream(fname):
        if eclipse != prev_eclipse:
            eclipses += 1
        if hitcount - last_flush_hitcount > 100000:
            flush_hits(hits, counts, False)
            last_flush_hitcount = hitcount
            progress(f"w{worker_id} [{fname.name}]: {eclipses} eclipses {hitcount} hits")

        for disk_index, point_index in pointset.query(disks, predicate="intersects").T:
            if has_nuv:
                hitcount += 1
                hits["n"][disk_index][point_index] += 1
            if has_fuv:
                hitcount += 1
                hits["f"][disk_index][point_index] += 1

    flush_hits(hits, counts, True)
    progress(f"w{worker_id} [{fname.name}]: {eclipses} eclipses {hitcount} hits, done.")


THE_POINTSET = None
THE_COUNTS = None
def crunch_disks_mp_wrapper(args):
    global THE_POINTSET
    global THE_COUNTS

    worker_id, fname = args
    crunch_disks(worker_id, fname, THE_POINTSET, THE_COUNTS)
    return worker_id


def process_sector(sector_fname, disks, output_writer):
    global THE_POINTSET
    global THE_COUNTS

    progress(f"loading {sector_fname}")
    with gzip.open(sector_fname, "rb") as fp:
        pointset = pickle.load(fp)

    progress("initializing counters")
    npoints = len(pointset.geometries)
    def z():
        return multiprocessing.Array("I", npoints)

    counts = {
        "n": [ z(), z(), z() ],
        "f": [ z(), z(), z() ],
    }

    THE_POINTSET = pointset
    THE_COUNTS = counts

    progress("counting intersections")
    with multiprocessing.Pool() as pool:
        for _ in pool.imap_unordered(
                crunch_disks_mp_wrapper, enumerate(disks), 1
        ):
            pass

    THE_POINTSET = None
    THE_COUNTS = None

    progress("constructing batch")
    # no need for locking at this point, we are single-threaded again
    batch = pa.RecordBatch.from_pydict({
        "gal_l":   [pt.x for pt in pointset.geometries],
        "gal_b":   [pt.y for pt in pointset.geometries],
        "nuv_350": counts["n"][0].get_obj(),
        "nuv_370": counts["n"][1].get_obj(),
        "nuv_400": counts["n"][2].get_obj(),
        "fuv_350": counts["f"][0].get_obj(),
        "fuv_370": counts["f"][1].get_obj(),
        "fuv_400": counts["f"][2].get_obj(),
    }, schema=RASTER_SCHEMA)

    progress("writing batch")
    output_writer.write_batch(batch)

    progress(f"finished with {sector_fname}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("disks_dir", type=Path)
    ap.add_argument("grid_dir", type=Path)
    ap.add_argument("output_parquet", type=Path)
    args = ap.parse_args()

    disks = list(args.disks_dir.glob("disks-*.gz"))

    with closing(parquet.ParquetWriter(
            args.output_parquet,
            RASTER_SCHEMA,
            use_dictionary = False,
            compression = "zstd",
            write_page_checksum=True
    )) as output_writer:
        for sector in args.grid_dir.glob("grid-tree-*.pickle.gz"):
            process_sector(sector, disks, output_writer)


main()
