"""
High-level stages of catalog generation
"""

from pathlib import Path
from typing import Optional

from glcat.cli import Band

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
    raise NotImplementedError


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
