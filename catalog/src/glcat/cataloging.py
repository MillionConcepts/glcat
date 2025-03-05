import csv
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd
from pyarrow import parquet

from gPhoton.aspect import aspect_tables
from gPhoton.eclipse import photomfile_path, expfile_path
from glcat.constants import Band, DEFAULT_APERTURES
from glcat.gphoton_ext import counts2mag, counts2flux


def accumulate_run_data(
    eclipse: int,
    leg: int,
    eclipse_dir: Path,
    aspect_dir: None | str,
    depth: int,
    aperture_sizes: Iterable[float],
):
    metadata = aspect_tables(eclipse, "metadata", aspect_dir=aspect_dir)[0]
    obstypes = set(s.as_py() for s in metadata["obstype"])
    obstype = obstypes.pop()
    assert not obstypes

    aperture_sizes = sorted(aperture_sizes)

    phot = {}
    expt = {}
    for band in [Band.NUV, Band.FUV]:
        expt[band] = pd.read_csv(
            eclipse_dir / expfile_path(
                eclipse,
                leg,
                band.name,
                mode = "direct",
                depth = depth,
                start = None,
            )
        )

        phot[band] = {}
        for mode in ["base", "forced"]:
            phot[band][mode] = {}
            for ap in aperture_sizes:
                phot[band][mode][ap] = parquet.read_table(
                    eclipse_dir / photomfile_path(
                        eclipse,
                        leg,
                        band.name,
                        mode = "direct",
                        depth = depth,
                        start = None,
                        aperture = ap,
                        suffix = mode,
                        ftype = "parquet"
                    )
                ).to_pandas()

    return {
        "eclipse": eclipse,
        "leg": leg,
        "obstype": obstype,
        "apertures": aperture_sizes,
        "expt": expt,
        "phot": phot,
    }



def make_catalog(
    data: dict[str, Any],
    band: Band,
    verbose: int = 0
) -> pd.DataFrame:
    """
    Produce a unified NUV/FUV photometric catalog of all sources
    detected in `band`, by combining the base photometry for that band
    with the *forced* photometry for the *other* band.
    """
    eclipse = data["eclipse"]
    obstype = data["obstype"]
    leg = data["leg"]
    apertures = data["apertures"]

    this_band = data["phot"][band]["base"]
    other_band = data["phot"][band.other]["forced"]

    # Index information is expected to be the same across all aperture
    # sizes, so we pull it from the photometry tables for an arbitrary
    # aperture size.
    ix_this_band = this_band[apertures[0]]
    ix_other_band = other_band[apertures[0]]

    nrows = len(ix_this_band)
    assert len(ix_other_band) == nrows

    # eclipse-wide metadata
    columns = [
        pd.Series(data=np.full(nrows, obstype), name="OBSTYPE"),
        pd.Series(data=np.full(nrows, eclipse), name="ECLIPSE"),
        pd.Series(data=np.full(nrows, leg), name="LEG"),
    ]
    columns.extend(
        pd.Series(name=f'APER_{i}', data=np.full(nrows, aper))
        for i, aper in enumerate(apertures)
    )
    # sky position index - same for both bands
    columns.append(ix_this_band[["ra", "dec"]])
    # extended source tag - only available for base photometry
    columns.append(ix_this_band["extended_source"].rename(f'{band}_EXTENDED'))

    # photometric centers for this band
    columns.append(ix_this_band[["xcenter", "ycenter"]].add_prefix(f"{band}_"))
    # photometric data for this band
    columns.extend(accumulate_photometry(data, band, "base", verbose))

    # photometric centers for the other band
    columns.append(
        ix_other_band[["xcenter", "ycenter"]].add_prefix(f"{band.other}_")
    )
    # photometric data for the other band
    columns.extend(accumulate_photometry(data, band.other, "forced", verbose))

    return pd.concat(columns, axis=1).rename(columns=str.upper)


def accumulate_photometry(
    data: dict[str, Any],
    band: Band,
    mode: str,
    verbose: int = 0
) -> pd.DataFrame:
    raw_expt = data["expt"][band]
    expt = raw_expt["expt"].sum()

    raw_phot = data["phot"][band][mode]
    nrows = len(next(iter(raw_phot.values())))

    cols = [pd.Series(data=np.full(nrows, expt), name=f'{band}_EXPT')]
    for i, aper in enumerate(data["apertures"]):
        cols.append(aper_photometry(expt, band, i, raw_phot[aper]))
    return cols


def aper_photometry(
    expt: float,
    band: Band,
    aper_ix: int,
    adat: pd.DataFrame,
) -> pd.DataFrame:
    cps = adat['aperture_sum'] / expt
    cps_err = np.sqrt(adat['aperture_sum']) / expt

    mag = counts2mag(cps, band)
    mag_err_upper = np.abs(counts2mag(cps - cps_err, band) - mag)
    mag_err_lower = np.abs(counts2mag(cps + cps_err, band) - mag)

    return pd.DataFrame({
        "SUM": adat["aperture_sum"],
        "EDGE": (adat["aperture_sum_edge"] != 0).astype(np.int8),
        "MASK": (adat["aperture_sum_mask"] != 0).astype(np.int8),
        "CPS": cps,
        "CPS_ERR": cps_err,
        "FLUX": counts2flux(cps, band),
        "FLUX_ERR": counts2flux(cps_err, band),
        "MAG": mag,
        "MAG_ERR_UPPER": mag_err_upper,
        "MAG_ERR_LOWER": mag_err_lower,
    }).add_prefix(f'{band}_').add_suffix(f'_A{aper_ix}')
