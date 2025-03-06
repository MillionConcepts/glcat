"""
High-level stages of catalog generation
"""

import shutil

from pathlib import Path
from typing import Optional, Sequence

from gPhoton.aspect import aspect_tables
from gPhoton.eclipse import raw6_path, photomfile_path, eclipse_prefix
from gPhoton.io.mast import retrieve_raw6
from gPhoton.pipeline import execute_pipeline
from gPhoton.reference import get_legs
from glcat.constants import Band, DEFAULT_APERTURES, DEFAULT_DEPTH
from glcat.cataloging import exposure_times_for_catalog, make_band_catalog


def download_raw(
    eclipse: int,
    *,
    eclipse_dir: Path,
    remote_eclipse_dir: Optional[Path] = None,
    bands: Sequence[Band] = [Band.NUV, Band.FUV],
    download: bool = True,
    recreate: bool = False,
    verbose: int = 0,
):
    for band in bands:
        raw6_relpath = raw6_path(eclipse, band.name, "direct")
        raw6_local = eclipse_dir / raw6_relpath
        raw6_local.parent.mkdir(parents=True, exist_ok=True)

        if raw6_local.exists():
            if verbose >= 2:
                print(f"eclipse {eclipse} band {band.name}:"
                      f" {raw6_local} already present")
            if not recreate:
                continue
        if remote_eclipse_dir is not None:
            raw6_remote = remote_eclipse_dir / raw6_relpath
            if raw6_remote.exists():
                if verbose >= 1:
                    print(f"eclipse {eclipse} band {band.name}:"
                          f" copying {raw6_local} from {raw6_remote}")
                shutil.copy(raw6_remote, raw6_local)
                continue
        if download:
            print(f"eclipse {eclipse} band {band.name}:"
                  f" downloading {raw6_local} from MAST")
            retrieve_raw6(eclipse, band.name, raw6_local)
        else:
            print(f"eclipse {eclipse} band {band.name}:"
                  f" {raw6_local} not available")


def base_photometry(
    eclipse: int,
    *,
    eclipse_dir: Path,
    remote_eclipse_dir: Optional[str] = None,
    aspect_dir: Optional[str] = None,
    bands: Sequence[Band] = [Band.NUV, Band.FUV],
    depth: Optional[int] = DEFAULT_DEPTH,
    aperture_sizes: list[float] = DEFAULT_APERTURES,
    recreate: bool = False,
    verbose: int = 0,
    parallel: int = 1,
    chunk_size: int = 1000000,
    share_memory: Optional[bool] = None,
):
    for band in bands:
        execute_pipeline(
            eclipse,
            band.name,
            threads = parallel,
            local_root = eclipse_dir,
            remote_root = remote_eclipse_dir,
            aspect_dir = aspect_dir,
            recreate = recreate,
            download = False,
            depth = depth,
            aperture_sizes = aperture_sizes,
            write = { "movie": True, "image": True },
            coregister_lightcurves = True,
            photometry_only = False,
            compression = "rice",
            suffix = "base",
            source_catalog_file = None,
            ftype = "parquet",
        )


def forced_photometry(
    eclipse: int,
    *,
    eclipse_dir: Path,
    remote_eclipse_dir: Optional[str] = None,
    aspect_dir: Optional[str] = None,
    bands: Sequence[Band] = [Band.NUV, Band.FUV],
    recreate: bool = False,
    depth: Optional[int] = DEFAULT_DEPTH,
    aperture_sizes: list[float] = DEFAULT_APERTURES,
    verbose: int = 0,
    parallel: int = 1,
    chunk_size: int = 1000000,
    share_memory: Optional[bool] = None,
):
    # The source catalog file is used only for source positions,
    # which will be the same for all apertures, so we arbitrarily
    # use the first one in the list.
    fp_src_aperture = aperture_sizes[0]

    for band in bands:
        for leg in get_legs(eclipse, aspect_dir=aspect_dir):
            fp_src = photomfile_path(
                eclipse,
                leg,
                band.other.name,
                "direct",
                depth = depth,
                start = None,
                aperture = fp_src_aperture,
                suffix = "base",
                ftype = "parquet",
            )
            execute_pipeline(
                eclipse,
                band.name,
                threads = parallel,
                local_root = eclipse_dir,
                remote_root = remote_eclipse_dir,
                aspect_dir = aspect_dir,
                recreate = recreate,
                download = False,
                depth = depth,
                aperture_sizes = aperture_sizes,
                write = { "movie": True, "image": True },
                coregister_lightcurves = True,
                photometry_only = True,
                compression = "rice",
                suffix = "forced",
                source_catalog_file = eclipse_dir / fp_src,
                ftype = "parquet",
            )


def band_catalog(
    eclipse: int,
    *,
    eclipse_dir: Path,
    aspect_dir: Optional[str] = None,
    bands: Sequence[Band] = [Band.NUV, Band.FUV],
    recreate: bool = False,
    depth: Optional[int] = DEFAULT_DEPTH,
    aperture_sizes: list[float] = DEFAULT_APERTURES,
    verbose: int = 0,
):
    metadata = aspect_tables(eclipse, "metadata", aspect_dir=aspect_dir)[0]
    obstypes = set(s.as_py() for s in metadata["obstype"])
    obstype = obstypes.pop()
    assert not obstypes

    for leg in get_legs(eclipse, aspect_dir=aspect_dir):
        # We always need the exposure time for both bands for each leg.
        exposure_times = exposure_times_for_catalog(
            eclipse, leg, depth=depth, eclipse_dir=eclipse_dir
        )

        for band in bands:
            d, p = eclipse_prefix(eclipse, band.name, "direct", False)
            catalog_path = eclipse_dir / d / f"{p}-{leg}-catalog.parquet"
            if recreate or not catalog_path.exists():
                make_band_catalog(
                    catalog_path,
                    eclipse = eclipse,
                    leg = leg,
                    band = band,
                    obstype = obstype,
                    eclipse_dir = eclipse_dir,
                    depth = depth,
                    aperture_sizes = aperture_sizes,
                    exposure_times = exposure_times,
                    verbose = verbose
                )


def merged_catalog(
    eclipse: int,
    *,
    eclipse_dir: Path,
    remote_eclipse_dir: Optional[str] = None,
    aspect_dir: Optional[str] = None,
    recreate: bool = False,
    verbose: int = 0,
):
    raise NotImplementedError
