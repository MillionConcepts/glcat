import csv
from pathlib import Path
from typing import Any

import numpy as np

# "import pyarrow as pa" does *not* make submodules available as pa.foo
import pyarrow
import pyarrow.compute
import pyarrow.csv
import pyarrow.parquet
pc = pyarrow.compute

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
        expt[band] = pc.sum(exp_data["expt"]).as_py()

    return expt


def photometry_for_catalog(
    eclipse: int,
    leg: int,
    depth: int,
    aperture_sizes: list[float],
    eclipse_dir: Path,
    band: Band,
    suffix: str,
) -> dict[float, pyarrow.Table]:
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
            ),
            columns = [
                "ra", "dec", "extended_source", "xcenter", "ycenter",
                "aperture_sum", "aperture_sum_mask", "aperture_sum_edge",
            ]
        ).sort_by([("ra","ascending"),
                   ("dec","ascending")])
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
    columns = {
        "OBSTYPE": np.full(nrows, obstype),
        "ECLIPSE": np.full(nrows, eclipse),
        "LEG":     np.full(nrows, leg),
    }
    for i, ap in enumerate(aperture_sizes):
        columns[f"APER_{i}"] = np.full(nrows, ap)

    # sky position of sources - same for both bands
    columns["RA"] = ix_base["ra"]
    columns["DEC"] = ix_base["dec"]
    # extended source tag - only available for base photometry
    columns[f"{b_base}_EXTENDED"] = ix_base["extended_source"]

    # photometric centers and exposure time for this band
    columns[f"{b_base}_XCENTER"] = ix_base["xcenter"]
    columns[f"{b_base}_YCENTER"] = ix_base["ycenter"]
    columns[f"{b_base}_EXPT"]    = np.full(nrows, exposure_times[b_base])

    # photometric data for this band
    for i, ap in enumerate(aperture_sizes):
        phot = ph_base[ap]
        assert pc.all(pc.equal(phot["ra"], ix_base["ra"]))
        assert pc.all(pc.equal(phot["dec"], ix_base["dec"]))
        columns.update(
            aper_photometry(exposure_times[b_base], b_base, i, phot)
        )

    # photometric centers for the other band
    columns[f"{b_forced}_XCENTER"] = ix_forced["xcenter"]
    columns[f"{b_forced}_YCENTER"] = ix_forced["ycenter"]
    columns[f"{b_forced}_EXPT"]    = np.full(nrows, exposure_times[b_forced])

    # photometric data for the other band
    for i, ap in enumerate(aperture_sizes):
        phot = ph_forced[ap]
        # intentionally comparing sky positions to the base index
        assert pc.all(pc.equal(phot["ra"], ix_base["ra"]))
        assert pc.all(pc.equal(phot["dec"], ix_base["dec"]))
        columns.update(
            aper_photometry(exposure_times[b_forced], b_forced, i, phot)
        )

    # and that's all! write out the table.
    catalog_table = pyarrow.Table.from_pydict(columns)
    pyarrow.parquet.write_table(
        catalog_table,
        catalog_path,
        # maximize interop with other parquet readers
        version="1.0",
        store_schema=False
    )


def aper_photometry(
    exposure_time: float,
    band: Band,
    aper_ix: int,
    phot: pyarrow.Table,
) -> dict[str, np.ndarray]:
    count = phot['aperture_sum'].to_numpy()
    cps = count / exposure_time
    cps_err = np.sqrt(count) / exposure_time

    mag = counts2mag(cps, band)
    mag_err_upper = np.abs(counts2mag(cps - cps_err, band) - mag)
    mag_err_lower = np.abs(counts2mag(cps + cps_err, band) - mag)

    return {
        f"{band}_SUM_A{aper_ix}":           count,
        f"{band}_EDGE_A{aper_ix}":          pc.not_equal(phot["aperture_sum_edge"], 0),
        f"{band}_MASK_A{aper_ix}":          pc.not_equal(phot["aperture_sum_mask"], 0),
        f"{band}_CPS_A{aper_ix}":           cps,
        f"{band}_CPS_ERR_A{aper_ix}":       cps_err,
        f"{band}_FLUX_A{aper_ix}":          counts2flux(cps, band),
        f"{band}_FLUX_ERR_A{aper_ix}":      counts2flux(cps_err, band),
        f"{band}_MAG_A{aper_ix}":           mag,
        f"{band}_MAG_ERR_UPPER_A{aper_ix}": mag_err_upper,
        f"{band}_MAG_ERR_LOWER_A{aper_ix}": mag_err_lower,
    }
