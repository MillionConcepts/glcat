import csv
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

import pyarrow
import pyarrow.csv
import pyarrow.parquet

from gPhoton.eclipse import photomfile_path, expfile_path
from glcat.constants import Band, DEFAULT_APERTURES
from glcat.gphoton_ext import counts2mag, counts2flux


def exposure_times_for_catalog(
    eclipse: int,
    leg: int,
    *,
    depth: int,
    eclipse_dir: Path,
) -> dict[Band, float]:
    """
    Read the exposure-time files for the current eclipse and leg.
    Compute the aggregate exposure times used by catalog generation.
    """

    # We only need the 'expt' column from the exposure-time files.
    options = pyarrow.csv.ConvertOptions(
        include_columns=["expt"],
        column_types={"expt": pyarrow.float64()}
    )

    # We always need the exposure time for both bands, whether or not
    # we're generating catalogs for both bands.
    expt = {}
    for band in [Band.NUV, Band.FUV]:
        exp_path = eclipse_dir / expfile_path(
            eclipse,
            leg,
            band.name,
            mode = "direct",
            depth = depth,
            start = None,
        )
        exp_data = pyarrow.csv.read_csv(exp_path, convert_options = options)
        expt[band] = exp_data["expt"].to_numpy().sum()

    return expt


def photometry_for_catalog(
    eclipse: int,
    leg: int,
    depth: int,
    aperture_sizes: list[float],
    eclipse_dir: Path,
    band: Band,
    suffix: str,
) -> dict[float, pd.DataFrame]:
    return {
        ap: pyarrow.parquet.read_table(
            eclipse_dir / photomfile_path(
                eclipse,
                leg,
                band.name,
                mode = "direct",
                depth = depth,
                start = None,
                aperture = ap,
                suffix = suffix,
                ftype = "parquet"
            )
        ).to_pandas()
        for ap in aperture_sizes
    }


def make_band_catalog(
    catalog_path: Path,
    *,
    eclipse: int,
    leg: int,
    band: Band,
    obstype: str,
    eclipse_dir: Path,
    depth: int,
    aperture_sizes: list[float],
    exposure_times: dict[Band, float],
    verbose: int
) -> None:
    """
    Produce a unified NUV/FUV photometric catalog of all sources
    detected in `band`, by combining the base photometry for that band
    with the *forced* photometry for the *other* band.

    Arguments:
      catalog_path   - Pathname where the catalog will be written.
      eclipse        - Eclipse number to be catalogued.
      leg            - Leg number to be catalogued.
      band           - Band to use for source positions.
      obstype        - Observation type for this eclipse.
      depth          - Number of seconds per movie frame.
      aperture_sizes - List of photometric aperture sizes to be included
                       in the catalog.  Assumed to be sorted.
      verbose        - Verbosity level.
    """

    b_base = band
    b_forced = band.other

    # Column prefixes
    p_base = f"{b_base}_"
    p_forced = f"{b_forced}_"

    # Photometry data for each band
    ph_base = photometry_for_catalog(
        eclipse, leg, depth, aperture_sizes, eclipse_dir, b_base, "base",
    )
    ph_forced = photometry_for_catalog(
        eclipse, leg, depth, aperture_sizes, eclipse_dir, b_forced, "forced",
    )

    # Index information is expected to be the same across all aperture
    # sizes, so we pull it from the photometry tables for an arbitrary
    # aperture size.
    ix_base = ph_base[aperture_sizes[0]]
    ix_forced = ph_forced[aperture_sizes[0]]

    nrows = len(ix_base)
    assert len(ix_forced) == nrows

    # construct the catalog table from:
    # eclipse-wide metadata
    columns = [
        pd.Series(data=np.full(nrows, obstype), name="OBSTYPE"),
        pd.Series(data=np.full(nrows, eclipse), name="ECLIPSE"),
        pd.Series(data=np.full(nrows, leg), name="LEG"),
    ]
    for i, ap in enumerate(aperture_sizes):
        columns.append(pd.Series(name=f'APER_{i}', data=np.full(nrows, ap)))

    # sky position of sources - same for both bands
    columns.append(ix_base[["ra", "dec"]])
    # extended source tag - only available for base photometry
    columns.append(ix_base["extended_source"].rename(f'{band}_EXTENDED'))

    # photometric centers and exposure time for this band
    columns.append(ix_base[["xcenter", "ycenter"]].add_prefix(p_base))
    columns.append(pd.Series(
        data=np.full(nrows, exposure_times[b_base]),
        name=f'{b_base}_EXPT'
    ))
    # photometric data for this band
    columns.extend(
        aper_photometry(exposure_times[b_base], b_base, i, ph_base[ap])
        for i, ap in enumerate(aperture_sizes)
    )

    # photometric centers for the other band
    columns.append(ix_forced[["xcenter", "ycenter"]].add_prefix(p_forced))
    columns.append(pd.Series(
        data=np.full(nrows, exposure_times[band.other]),
        name=f'{band.other}_EXPT'
    ))

    # photometric data for the other band
    columns.extend(
        aper_photometry(exposure_times[b_forced], b_forced, i, ph_forced[ap])
        for i, ap in enumerate(aperture_sizes)
    )

    # and that's all
    pd.concat(columns, axis=1).rename(columns=str.upper).to_parquet(
        catalog_path
    )


def aper_photometry(
    exposure_time: float,
    band: Band,
    aper_ix: int,
    phot: pd.DataFrame,
) -> pd.DataFrame:
    cps = phot['aperture_sum'] / exposure_time
    cps_err = np.sqrt(phot['aperture_sum']) / exposure_time

    mag = counts2mag(cps, band)
    mag_err_upper = np.abs(counts2mag(cps - cps_err, band) - mag)
    mag_err_lower = np.abs(counts2mag(cps + cps_err, band) - mag)

    return pd.DataFrame({
        "SUM": phot["aperture_sum"],
        "EDGE": (phot["aperture_sum_edge"] != 0).astype(np.int8),
        "MASK": (phot["aperture_sum_mask"] != 0).astype(np.int8),
        "CPS": cps,
        "CPS_ERR": cps_err,
        "FLUX": counts2flux(cps, band),
        "FLUX_ERR": counts2flux(cps_err, band),
        "MAG": mag,
        "MAG_ERR_UPPER": mag_err_upper,
        "MAG_ERR_LOWER": mag_err_lower,
    }).add_prefix(f'{band}_').add_suffix(f'_A{aper_ix}')
