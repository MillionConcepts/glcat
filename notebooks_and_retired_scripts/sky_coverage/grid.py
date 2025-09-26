"""Phase 2 of sky coverage map generation.  Generate a set of
spatially indexed point sets (as shapely STRtree objects) which
collectively form a (nearly) equally spaced grid on the celestial
sphere, with average separation between nearest neighbors of 1 minute
of arc.

Each of the point sets covers an octant of the sky, and is written out
to a separate file as a pickle.  We do it this way because the complete
grid on the celestial sphere is too large to hold in memory all at once.

This phase works strictly with abstract coordinates; each grid point
is (lon, lat) where -180° < lon < 180° and -90° < lat < 90° and the
origin is *unspecified.*
"""

import argparse
import gzip
import os
import pickle
import pickletools
import sys
import time

from datetime import timedelta
from math import pi as PI
from pathlib import Path

import numpy as np
import shapely

# Conversions between radians and degrees
DEG_TO_RAD = PI/180.0
RAD_TO_DEG = 180.0/PI

# Target angular resolution for the map is 1′
ANGULAR_RES = 1/60


START_TIME = None
def progress(msg: str) -> None:
    now = time.monotonic()

    global START_TIME
    if START_TIME is None:
        START_TIME = now

    elapsed = timedelta(seconds=now - START_TIME)
    sys.stderr.write(f"[{elapsed}]  {msg}\n")


def distribute_points_on_sphere(*, n=None, k=None):
    """
    Produce points uniformly distributed over the surface of a sphere.
    "Uniformly distributed" means the average distance between each
    point and its nearest neighbor is as close to constant as possible.

    Caller must supply either 'n', the number of points to generate,
    or 'k', the desired average distance between nearest neighbors in
    degrees.

    The algorithm is taken from
    https://extremelearning.com.au/how-to-evenly-distribute-points-on-a-sphere-more-effectively-than-the-canonical-fibonacci-lattice/
    (see section 3 re optimizing the _average_ rather than _minimum_
    distance between points).
    """

    match (n, k):
        case (_, None):
            pass
        case (None, _):
            # As we are distributing points over a two-dimensional surface,
            # the average distance between nearest neighbors is inversely
            # proportional to the square root of the number of points.
            # Empirically, k = 196.62928770516 / √n -- there is probably
            # a nice closed-form expression for the constant but I don't
            # rememeber enough trig to find it.  Solving for n gives:
            n = 38663.0767834385 / (k*k)
        case (_, _):
            raise ValueError("must supply one of n or k, but not both")

    # The algorithm below breaks down if n < 2.  In these cases, pick
    # points that avoid the poles and the ±180° coordinate discontinuity.
    if n == 0:
        return np.array([]), np.array([])
    if n == 1:
        return np.array([0.0]), np.array([0.0])
    if n == 2:
        return np.array([-90, 90]), np.array([0, 0])

    epsilon = 0.36

    i = np.arange(n)
    lon = (4 * 180 / (1 + np.sqrt(5))) * i % 360 - 180
    lat = np.arccos(np.minimum(1.0, np.maximum(-1.0,
        1 - 2*(i + epsilon)/(n - 1 + 2*epsilon)
    ))) * RAD_TO_DEG - 90

    return lon, lat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("outdir", type=Path)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    progress("computing points")
    lon, lat = distribute_points_on_sphere(k=ANGULAR_RES)

    # At an angular resolution of 1′ the point grid is too large for a
    # workstation with 32G of RAM to process efficiently (it only goes
    # into swap a _little_, but that's still the difference between
    # taking 15 seconds and taking _seven minutes_ to build an R-tree).
    # Divide it into eight chunks and write each chunk to disk.
    # The next stage will process each chunk separately.
    chunk = 0
    for lat_min, lat_max in [(-90, 0), (0, 91)]:
        progress(f"lat [{lat_min}, {lat_max})...")
        lat_range = (lat >= lat_min) & (lat < lat_max)
        for lon_min, lon_max in [(-180, -90), (-90, 0), (0, 90), (90, 181)]:
            progress(f"lon [{lon_min}, {lon_max}): building tree...")
            lon_range = (lon >= lon_min) & (lon < lon_max)
            tree = shapely.STRtree(shapely.points(
                lon[lat_range & lon_range],
                lat[lat_range & lon_range]
            ))

            with gzip.open(
                    args.outdir / f"grid-tree-{chunk}.pickle.gz", "wb"
            ) as fp:
                progress(f"lon [{lon_min}, {lon_max}): pickling tree...")
                pkl = pickle.dumps(tree, -1)
                tree = None
                progress(f"lon [{lon_min}, {lon_max}): squeezing pickle...")
                pkl = pickletools.optimize(pkl)
                progress(f"lon [{lon_min}, {lon_max}): zipping pickle...")
                fp.write(pkl)
                pkl = None
                progress(f"lon [{lon_min}, {lon_max}): writing...")
                fp.flush()
                os.fsync(fp.fileno())

            chunk += 1


main()
