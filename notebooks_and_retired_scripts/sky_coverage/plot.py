"""Phase 4 of sky coverage map generation:  Turn the "raster" generated
by phase 3 into human-meaningful maps."""

import argparse

from pathlib import Path
from math import log10, pi as PI

import colorcet as cc
import numpy as np
import matplotlib.pyplot as plt

from cartopy import crs
from matplotlib.ticker import FuncFormatter
from pyarrow import parquet


CEL_SPHERE = crs.Globe(
    ellipse=None,
    # these axes produce a sphere whose projections have 1 map unit = 1 degree,
    # which is what we need
    semimajor_axis=180/PI,
    semiminor_axis=180/PI,
)
CRS_GALACTIC = crs.Geodetic(globe=CEL_SPHERE)
PC_GALACTIC = crs.PlateCarree(globe=CEL_SPHERE)

# roughly evenly spaced breaks between 1s and 20h
# see scale() re why log10(1.5) is labeled "1 s"
EXPOSURE_LBLS = {
    log10(1.5):   "1 s",
    log10(5):     "5 s",
    log10(20):    "20 s",
    log10(60):    "1 m",
    log10(300):   "5 m",
    log10(1200):  "20 m",
    log10(3600):  "1 h",
    log10(18000): "5 h",
    log10(72000): "20 h",
}

EXPOSURE_TICKS = sorted(EXPOSURE_LBLS.keys())
EXPOSURE_FMT = FuncFormatter(lambda x, _: EXPOSURE_LBLS[x])


def scale(exposure):
    """
    put EXPOSURE on a log scale, except that 0 is left at 0 and 1
    is mapped to log10(1.5).  Values >= 2 are simply run through
    log10().  (N.B. It is known that all values of EXPOSURE are
    integers >= 0.)
    """
    normal_ = exposure > 1
    one_ = exposure == 1

    scaled = np.zeros_like(exposure)
    scaled[one_] = log10(1.5)
    scaled[normal_] = np.log10(exposure[normal_])
    return scaled


def plot_exposure(fig, ax, raster, col):
    l = raster["gal_l"]
    b = raster["gal_b"]
    t = scale(raster[col])
    pos = t > 0

    ax.set_title(col)
    return ax.scatter(
        l[pos], b[pos], c=t[pos],
        marker='o',
        lw=0,
        s=(72./fig.dpi)**2,
        transform=PC_GALACTIC,
        cmap=cc.m_CET_L6_r,
    )


def plot_all(base_path, raster, width, height):
    for col in raster.columns:
        if col.startswith("gal_"):
            continue

        path = base_path.with_stem(base_path.stem + "_" + col)
        fig = plt.figure(layout="constrained", figsize=(width, height))
        try:
            ax = fig.subplots(1, 1, subplot_kw={
                "xlim": (-180, 180),
                "ylim": (-90, 90),
                "projection": crs.Mollweide(globe=CEL_SPHERE)
            })
            im = plot_exposure(fig, ax, raster, col)
            bar = fig.colorbar(
                im,
                shrink=0.6,
                ax=ax,
                ticks=EXPOSURE_TICKS,
                format=EXPOSURE_FMT,
                label="Exposure time",
            )
            fig.savefig(path)
        finally:
            plt.close(fig)



def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("plot_base", type=Path)
    ap.add_argument("raster", type=Path)
    ap.add_argument("width", type=float)
    ap.add_argument("height", type=float)
    a = ap.parse_args()

    raster = parquet.read_table(a.raster).to_pandas()
    plot_all(a.plot_base, raster, a.width, a.height)


main()
