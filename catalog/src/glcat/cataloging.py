import csv
from pathlib import Path
from typing import Any

import numpy as np
from datetime import datetime, timezone

# "import pyarrow as pa" does *not* make submodules available as pa.foo
import pyarrow
import pyarrow.compute
import pyarrow.csv
import pyarrow.parquet
pc = pyarrow.compute

from gPhoton.eclipse import photomfile_path, expfile_path
from glcat.constants import Band, DEFAULT_APERTURES
from glcat.gphoton_ext import counts2mag, counts2flux
from glcat.modelphotometry import model_photometry


def exposure_times_for_catalog(
    eclipse: int,
    leg: int,
    *,
    depth: int,
    eclipse_dir: Path,
    aperture_size: float,
) -> dict[Band, float]:
    """
    Read the exposure-time metadata for the current eclipse and leg
    from the photometry tables. This *used* to be a standalone csv.
    Compute the aggregate exposure times used by catalog generation.
    Default to photom table for first aperture in aperture list. """
    expt = {}
    for band in [Band.NUV, Band.FUV]:
        photom_path = eclipse_dir / photomfile_path(
            eclipse,
            leg,
            band.name,
            mode = "direct",
            depth = depth,
            start = None,
            aperture = aperture_size,
            ftype = "parquet",
            suffix = "base",
        )
        try:
            photom_file = pyarrow.parquet.ParquetFile(photom_path)
            metadata = photom_file.schema_arrow.metadata
            if metadata is None or b"XPOSURE" not in metadata:
                raise ValueError(f"'xposure' not in metadata of {photom_path}")
            xposure_str = float(metadata[b"XPOSURE"].decode("utf-8"))
            expt[band] = float(xposure_str)
        except FileNotFoundError:
            raise FileNotFoundError(
                f"eclipse {eclipse} band {band}: photom table {photom_path} not found"
            )

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
    photometry: dict[float, pyarrow.Table] = {}
    for ap in aperture_sizes:
        src = eclipse_dir / photomfile_path(
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
        columns = [
            "ra", "dec", "extended_source", "xcenter", "ycenter",
            "aperture_sum", "artifact_flag"
        ]
        if suffix == "base":
            columns.extend(["ya_aperture_sum", "stdcolrow_aperture_sum","q_aperture_sum"])

        if src.exists():
            photometry[ap] = pyarrow.parquet.read_table(
                src,
                columns=columns
            ).sort_by([("ra", "ascending"),
                       ("dec", "ascending")])
        else:
            print(f"eclipse {eclipse} band {band} aperture {ap}: "
                  f"photomfile {src} not found")
            photometry[ap] = pyarrow.Table.from_pydict({
                col: [] for col in columns
            })
    return photometry


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
    with the *forced* photometry for the *other* band. Also run model
    based photometry.

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
    if nrows == 0:
        return
    assert len(ix_forced) in (0, nrows), f"for {leg}, {len(ix_forced)} vs {len(ix_base)}"

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
        columns.update(
            aper_alt(b_base, i, phot)
        )

    # add model photometry columns for both base and forced bands
    columns.update(model_photometry(aperture_sizes, columns, b_base))

    if len(ix_forced) == nrows:
        # photometric centers for the other band
        columns[f"{b_forced}_XCENTER"] = ix_forced["xcenter"]
        columns[f"{b_forced}_YCENTER"] = ix_forced["ycenter"]
        columns[f"{b_forced}_EXPT"] = np.full(nrows, exposure_times[b_forced])

        # photometric data for the other band
        for i, ap in enumerate(aperture_sizes):
            phot = ph_forced[ap]
            # intentionally comparing sky positions to the base index
            assert pc.all(pc.equal(phot["ra"], ix_base["ra"]))
            assert pc.all(pc.equal(phot["dec"], ix_base["dec"]))
            columns.update(
                aper_photometry(exposure_times[b_forced], b_forced, i, phot)
            )

        # add model photometry columns for forced band
        columns.update(model_photometry(aperture_sizes, columns, b_forced))

    else:
        # no data available for the other band, fill in blank columns
        nulls = pyarrow.nulls
        columns[f"{b_forced}_XCENTER"] = nulls(nrows)
        columns[f"{b_forced}_YCENTER"] = nulls(nrows)
        columns[f"{b_forced}_EXPT"]    = nulls(nrows)
        for i, ap in enumerate(aperture_sizes):
            columns[f"{b_forced}_SUM_A{i}"]           = nulls(nrows)
            columns[f"{b_forced}_FLAG_A{i}"]          = nulls(nrows)
            columns[f"{b_forced}_CPS_A{i}"]           = nulls(nrows)
            columns[f"{b_forced}_CPS_ERR_A{i}"]       = nulls(nrows)
            columns[f"{b_forced}_FLUX_A{i}"]          = nulls(nrows)
            columns[f"{b_forced}_FLUX_ERR_A{i}"]      = nulls(nrows)
            columns[f"{b_forced}_MAG_A{i}"]           = nulls(nrows)
            columns[f"{b_forced}_MAG_ERR_UPPER_A{i}"] = nulls(nrows)
            columns[f"{b_forced}_MAG_ERR_LOWER_A{i}"] = nulls(nrows)
            # model photometry cols
            columns[f"{b_forced}_MDL_SE_A{i}"]        = nulls(nrows)
            columns[f"{b_forced}_MDL_RSL_A{i}"]       = nulls(nrows)
            for name in ["CPS", "SIGMA", "BKG_CPS", "CPS_SE", "SIGMA_SE", "BKG_CPS_SE", "DET"]:
                columns[f"{b_forced}_MDL_{name}"]     = nulls(nrows)

    # add a check that there's no duplicate col names, to prevent mayhem
    if len(set(columns.keys())) != len(columns):
        raise ValueError("duplicate column names detected after model photometry")

    # index column
    row_index = np.arange(nrows)
    index = np.array([
        f"{eclipse:05}_{leg:02}_{idx:06}_{band}"
        for idx in row_index])
    columns["BAND_INDEX"] = index

    # and that's all! write out the table.
    catalog_table = pyarrow.Table.from_pydict(columns)

    metadata = {
        b"TELESCOP": b"GALEX",
        b"ECLIPSE": str(eclipse).encode(),
        b"LEG": str(leg).encode(),
        b"BANDNAME": str(band).encode(),
        b"BAND": str(1 if band == "NUV" else 2).encode(),
        b"ORIGIN": b"Million Concepts",
        b"DATE": datetime.now(timezone.utc).replace(microsecond=0).isoformat().encode(),
        b"TIMESYS": b"UTC",
        b"VERSION": f"GLCAT_1.0",
    }
    catalog_table = catalog_table.replace_schema_metadata(metadata)

    pyarrow.parquet.write_table(
        catalog_table,
        catalog_path,
        # maximize interop with other parquet readers
        version="1.0",
        store_schema=True
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

    # dividing artifact flags into individual flags
    artifact_flags = phot['artifact_flag'].to_numpy()
    hotspot_flag = ((artifact_flags & 1) > 0).astype(int)
    ghost_flag = ((artifact_flags & 2) > 0).astype(int)
    soft_edge_flag = ((artifact_flags & 4) > 0).astype(int)
    hard_edge_flag = ((artifact_flags & 8) > 0).astype(int)

    return {
        f"{band}_SUM_A{aper_ix}": count,
        f"{band}_CPS_A{aper_ix}": cps,
        f"{band}_CPS_ERR_A{aper_ix}": cps_err,
        f"{band}_FLUX_A{aper_ix}": counts2flux(cps, band),
        f"{band}_FLUX_ERR_A{aper_ix}": counts2flux(cps_err, band),
        f"{band}_MAG_A{aper_ix}": mag,
        f"{band}_MAG_ERR_UPPER_A{aper_ix}": mag_err_upper,
        f"{band}_MAG_ERR_LOWER_A{aper_ix}": mag_err_lower,
        f"{band}_HOTSPOT_FLAG_A{aper_ix}": hotspot_flag,
        f"{band}_GHOST_FLAG_A{aper_ix}": ghost_flag,
        f"{band}_SOFTEDGE_FLAG_A{aper_ix}": soft_edge_flag,
        f"{band}_HARDEDGE_FLAG_A{aper_ix}": hard_edge_flag,

    }


def aper_alt(
    band: Band,
    aper_ix: int,
    phot: pyarrow.Table,
) -> dict[str, np.ndarray]:

    ya_photom = phot['ya_aperture_sum'].to_numpy()
    stdcolrow_photom = phot['stdcolrow_aperture_sum'].to_numpy()
    q_photom = phot['q_aperture_sum'].to_numpy()

    return {
        f"{band}_YA_A{aper_ix}": ya_photom,
        f"{band}_SIGDISP_A{aper_ix}": stdcolrow_photom,
        f"{band}_Q_A{aper_ix}": q_photom
    }

