"""
Tool for generating star catalogues from the raw data collected by the
GALEX ultraviolet telescope.

For details of each available action use 'glcat <action> --help'.
"""

import argparse
import enum
import sys

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, NoReturn


class Band(enum.Flag):
    """Indicates which band(s) of the raw data to process in a particular
       operation."""
    NUV = 0x01
    FUV = 0x02
    ALL = 0x03

    @classmethod
    def parse(cls, s: str) -> 'Band':
        """Parse a string as a Band value or comma-separated sequence
           of Band values which are or-ed together.  Parsing is case
           insensitive and accepts several convenient aliases."""

        val = cls(0)
        for token in s.lower().split(','):
            token = token.strip()
            if token in ("nuv", "near"):
                val |= cls.NUV
            elif token in ("fuv", "far"):
                val |= cls.FUV
            elif token in ("both", "all"):
                val |= cls.ALL
            else:
                raise ValueError(f"unrecognized frequency band {token!r}")
        return val


class Stage(enum.Flag):
    """Indicates which stages of processing to attempt in this run."""
    DOWNLOAD          = 0x01
    BASE_PHOTOMETRY   = 0x02
    FORCED_PHOTOMETRY = 0x04
    BAND_CATALOG      = 0x08
    MERGED_CATALOG    = 0x10
    ALL               = 0x1F

    @classmethod
    def parse(cls, s: str) -> 'Stage':
        """Parse a string as a Stage value or comma-separated sequence
           of Stage values which are or-ed together.  Parsing is case
           insensitive and accepts several convenient aliases."""

        val = cls(0)
        for token in s.lower().split(','):
            token = token.strip()
            if token == "all":
                val |= cls.ALL
            elif token in ("download", "down", "dl"):
                val |= cls.DOWNLOAD
            elif token in ("base_photometry", "base_photo", "base"):
                val |= cls.BASE_PHOTOMETRY
            elif token in ("forced_photometry", "force_photometry",
                           "forced_photo", "force_photo",
                           "forced", "force"):
                val |= cls.FORCED_PHOTOMETRY
            elif token in ("band_catalog", "band_cat", "band"):
                val |= cls.BAND_CATALOG
            elif token in ("merged_catalog", "merge_catalog",
                           "merged", "merge"):
                val |= cls.MERGED_CATALOG
            else:
                raise ValueError(f"unrecognized stage {stage!r}")
        return val


@dataclass
class CLIOptions:
    """Command line options after resolving all defaults and
       converting from strings.
    """
    eclipses: list[int]
    bands: Band
    stages: Stage
    local_root: Path
    remote_root: Optional[str]
    download: bool
    recreate: bool
    verbose: int
    parallel: int
    chunk_size: int
    share_memory: Optional[bool]
    aspect_dir: Path

    @classmethod
    def from_args(cls, args: Optional[Iterable[str]] = None) -> 'CLIOptions':
        ap = argparse.ArgumentParser(description=__doc__)
        # note: nothing gets a short option yet (except --verbose and
        # --parallel, which have standard short options) because
        # I need to see the complete set of options across *all*
        # subcommands to know what *deserves* a short option
        ap.add_argument(
            # https://github.com/python/cpython/issues/59330 and its
            # many duplicates
            dest="local_root",
            metavar="local-catalog-root",
            type=Path,
            help="Root of a local directory tree containing catalog files,"
            " cached raw6 and aspect files, etc.  Must be writable."
        )
        ap.add_argument(
            "eclipse", type=int, nargs="+",
            help="Eclipse numbers of eclipses to process.",
        )
        ap.add_argument(
            "--band", type=Band.parse, default=Band.ALL, metavar="NUV|FUV|both",
            help="comma-separated list of spectral bands to process",
        )
        ap.add_argument(
            "--stage", type=Stage.parse, default=Stage.ALL,
            help="comma-separated list of processing stages (default: all)"
        )
        ap.add_argument(
            "--remote-root", metavar="DIR",
            help="alternate location to check for existing raw6 and photonlist files"
        )
        ap.add_argument(
            "--download", action="store_true",
            help="download raw6 files from MAST as necessary",
        )
        ap.add_argument(
            "--recreate", action="store_true",
            help="erase and re-create any outputs that already exist",
        )
        ap.add_argument(
            "--verbose", "-v", action="count",
            help="be more verbose (repeat for even more detail)",
        )
        ap.add_argument(
            "--parallel", "-j", type=int, default=1, metavar="TASKS",
            help="maximum number of parallel tasks to run"
        )
        ap.add_argument(
            "--chunk-size", type=int, default=1000000, metavar="PHOTONS",
            help="maximum number of photons per photonpipe chunk"
            " (default: 1,000,000)"
        )
        ap.add_argument(
            "--aspect-dir", type=Path, default=None,
            help="find aspect files in this directory instead of"
            " <local_root>/aspect"
        )
        ap.add_argument(
            "--force-enable-shared-memory",
            action="store_const", dest="share_memory", const=True, default=None,
            help="force use of shared memory even when running single-threaded"
            " (mostly useful for debugging)"
        )
        ap.add_argument(
            "--force-disable-shared-memory",
            action="store_const", dest="share_memory", const=False, default=None,
            help="force non-use of shared memory even when running multithreaded"
            " (mostly useful for debugging)"
        )

        ns = ap.parse_args(args)

        if 'aspect_dir' in ns and ns.aspect_dir is not None:
            aspect_dir = ns.aspect_dir
        else:
            aspect_dir = ns.local_root / "aspect"

        return cls(
            eclipses = ns.eclipse,
            stages = ns.stage,
            bands = ns.band,
            local_root = ns.local_root,
            remote_root = ns.remote_root,
            download = ns.download,
            recreate = ns.recreate,
            verbose = ns.verbose,
            parallel = ns.parallel,
            chunk_size = ns.chunk_size,
            share_memory = ns.share_memory,
            aspect_dir = aspect_dir,
        )


def process_eclipse(eclipse: int, options: CLIOptions):
    # This import has to be at the function level to avoid a circular
    # dependency.
    from glcat import stages

    if options.stages & Stage.DOWNLOAD:
        stages.download_raw(
            eclipse,
            local_root=options.local_root,
            remote_root=options.remote_root,
            bands=options.bands,
            download=options.download,
            recreate=options.recreate,
            verbose=options.verbose,
        )
    if options.stages & Stage.BASE_PHOTOMETRY:
        stages.base_photometry(
            eclipse,
            local_root=options.local_root,
            remote_root=options.remote_root,
            aspect_dir=options.aspect_dir,
            bands=options.bands,
            recreate=options.recreate,
            verbose=options.verbose,
            parallel=options.parallel,
            chunk_size=options.chunk_size,
            share_memory=options.share_memory
        )
    if options.stages & Stage.FORCED_PHOTOMETRY:
        stages.forced_photometry(
            eclipse,
            local_root=options.local_root,
            remote_root=options.remote_root,
            aspect_dir=options.aspect_dir,
            bands=options.bands,
            recreate=options.recreate,
            verbose=options.verbose,
            parallel=options.parallel,
            chunk_size=options.chunk_size,
            share_memory=options.share_memory
        )
    if options.stages & Stage.BAND_CATALOG:
        stages.band_catalog(
            eclipse,
            local_root=options.local_root,
            remote_root=options.remote_root,
            aspect_dir=options.aspect_dir,
            bands=options.bands,
            recreate=options.recreate,
            verbose=options.verbose,
        )
    if options.stages & Stage.MERGED_CATALOG:
        stages.merged_catalog(
            eclipse,
            local_root=options.local_root,
            remote_root=options.remote_root,
            aspect_dir=options.aspect_dir,
            recreate=options.recreate,
            verbose=options.verbose,
        )


def main() -> NoReturn:
    options = CLIOptions.from_args()

    # TODO: cross-stage/cross-eclipse parallelism?
    # spin up a multiprocessing Pool and farm out invocations
    # of process_eclipse.  requires figuring out how to split
    # options.parallel into eclipse-level and chunk-level parallelism.
    for eclipse in options.eclipses:
        process_eclipse(eclipse, options)


if __name__ == "__main__":
    main()
