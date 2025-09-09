import argparse
import os

import matplotlib.pyplot as plt
import shapely

from cartopy import crs
from pyarrow import parquet
from math import pi as PI

CEL_SPHERE = crs.Globe(
    ellipse=None,
    # these axes produce a sphere whose projections have 1 map unit = 1 degree,
    # which is what aperture_disks etc need
    semimajor_axis=180/PI,
    semiminor_axis=180/PI,
)
PC_EQUATORIAL = crs.PlateCarree(globe=CEL_SPHERE)

COLORS = [
    '#1b9e77',
    '#d95f02',
    '#7570b3',
    '#e7298a',
    '#66a61e',
    '#e6ab02',
    '#a6761d',
    '#666666',
]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("table")
    ap.add_argument("outdir")
    args = ap.parse_args()

    md = parquet.read_table(args.table).to_pandas()
    os.makedirs(args.outdir, exist_ok=True)

    fig = None
    ax = None
    n = 0
    figs = 0
    for row in md.itertuples(index=False):
        if (row.ra_min != row.ra_min
            or row.ra_max != row.ra_max
            or row.dec_min != row.dec_min
            or row.dec_max != row.dec_max):
            continue

        if n == 0:
            if fig is not None:
                figs += 1
                fig.savefig(f"{args.outdir}/bbox-sanity-{figs:03}.png")
                plt.close(fig)
            fig = plt.figure()
            ax = fig.add_subplot(projection=PC_EQUATORIAL)
            ax.gridlines()
            ax.set_global()

        n = (n + 1) % 100

        ax.add_geometries(
            [shapely.box(min(row.ra_min, row.ra_max),
                         row.dec_min,
                         max(row.ra_min, row.ra_max),
                         row.dec_max)],
            facecolor='none',
            edgecolor=COLORS[n % len(COLORS)],
            crs=PC_EQUATORIAL,
        )

    if fig is not None:
        figs += 1
        fig.savefig(f"{args.outdir}/bbox-sanity-{figs:03}.png")
        plt.close(fig)

main()
