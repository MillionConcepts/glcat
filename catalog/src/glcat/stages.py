"""
High-level stages of catalog generation
"""

from pathlib import Path
from typing import Optional

from gPhoton.pipeline import execute_pipeline
from glcat.constants import Band, DEFAULT_APERTURES, DEFAULT_DEPTH


def download_raw(
    eclipse: int,
    *,
    local_root: Path,
    remote_root: Optional[str] = None,
    bands: Band = Band.ALL,
    download: bool = True,
    recreate: bool = False,
    verbose: int = 0,
):
    raise NotImplementedError


def base_photometry(
    eclipse: int,
    *,
    local_root: Path,
    remote_root: Optional[str] = None,
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
            local_root = local_root,
            remote_root = remote_root,
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
    local_root: Path,
    remote_root: Optional[str] = None,
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
    local_root: Path,
    remote_root: Optional[str] = None,
    aspect_dir: Optional[str] = None,
    bands: Band = Band.ALL,
    recreate: bool = False,
    verbose: int = 0,
):
    raise NotImplementedError


def merged_catalog(
    eclipse: int,
    *,
    local_root: Path,
    remote_root: Optional[str] = None,
    aspect_dir: Optional[str] = None,
    recreate: bool = False,
    verbose: int = 0,
):
    raise NotImplementedError
