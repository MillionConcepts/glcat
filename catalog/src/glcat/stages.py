"""
High-level stages of catalog generation
"""

from pathlib import Path
from typing import Optional

from gPhoton.pipeline import execute_pipeline
from glcat.cli import Band

# TODO Should these be command line parameters or pulled from eclipse
# metadata or something?
APERTURE_SIZES = [ 1.5, 2.3, 3.8, 6.0, 9.0, 12.8, 17.3, 30., 60., 90., ]
DEPTH = 120

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
            depth = DEPTH,
            aperture_sizes = APERTURE_SIZES,
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
