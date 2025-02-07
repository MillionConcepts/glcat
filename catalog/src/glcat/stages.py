"""
High-level stages of catalog generation
"""

import shutil

from pathlib import Path
from typing import Optional

from gPhoton.io.mast import retrieve_raw6
from gPhoton.pipeline import execute_pipeline
from gPhoton.eclipse import raw6_path
from glcat.constants import Band, DEFAULT_APERTURES, DEFAULT_DEPTH


def download_raw(
    eclipse: int,
    *,
    eclipse_dir: Path,
    remote_eclipse_dir: Optional[Path] = None,
    bands: Band = Band.ALL,
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
    bands: Band = Band.ALL,
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
            suffix = band.suffix,
            source_catalog_file = None
        )


def forced_photometry(
    eclipse: int,
    *,
    eclipse_dir: Path,
    remote_eclipse_dir: Optional[str] = None,
    aspect_dir: Optional[str] = None,
    bands: Band = Band.ALL,
    recreate: bool = False,
    depth: Optional[int] = DEFAULT_DEPTH,
    aperture_sizes: list[float] = DEFAULT_APERTURES,
    verbose: int = 0,
    parallel: int = 1,
    chunk_size: int = 1000000,
    share_memory: Optional[bool] = None,
):
    raise NotImplementedError


def band_catalog(
    eclipse: int,
    *,
    eclipse_dir: Path,
    remote_eclipse_dir: Optional[str] = None,
    aspect_dir: Optional[str] = None,
    bands: Band = Band.ALL,
    recreate: bool = False,
    verbose: int = 0,
):
    raise NotImplementedError


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
