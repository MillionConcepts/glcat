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

from gPhoton.reference import check_eclipse
from glcat import stages
from glcat.constants import (
    Band,
    parse_bands,
    ALL_APERTURES,
    DEFAULT_APERTURES,
    DEFAULT_DEPTH,
)


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
    depth: Optional[int]
    apertures: list[float]
    download: bool
    recreate: bool
    verbose: int
    parallel: int
    chunk_size: int
    share_memory: Optional[bool]
    catalog_root: Path
    aspect_dir: Path
    eclipse_dir: Path
    remote_eclipse_dir: Optional[Path]
    print_tracebacks: bool

    @classmethod
    def from_args(cls, args: Optional[Iterable[str]] = None) -> 'CLIOptions':

        def parse_depth(arg: str) -> Optional[int]:
            arg = arg.lower()
            if arg in ("-", "inf", "infinity", "unlimited"):
                return None
            return int(arg)

        def parse_apertures(arg: str) -> list[float]:
            arg = arg.lower()
            if arg in ("all", "mission"):
                return ALL_APERTURES
            if arg == "default":
                return DEFAULT_APERTURES
            return sorted(set(float(a.strip()) for a in arg.split(",")))

        ap = argparse.ArgumentParser(description=__doc__)
        # note: nothing gets a short option yet (except --verbose and
        # --parallel, which have standard short options) because
        # I need to see the complete set of options across *all*
        # subcommands to know what *deserves* a short option
        ap.add_argument(
            # https://github.com/python/cpython/issues/59330 and its
            # many duplicates
            dest="catalog_root",
            metavar="catalog-root",
            type=Path,
            help="Root of a local directory tree containing catalog files,"
            " cached raw6 and aspect files, etc.  Must be writable."
        )
        ap.add_argument(
            "eclipse", type=int, nargs="+",
            help="Eclipse numbers of eclipses to process.",
        )
        ap.add_argument(
            "--band", type=parse_bands, default=[Band.NUV, Band.FUV],
            metavar="NUV|FUV|both",
            help="comma-separated list of spectral bands to process",
        )
        ap.add_argument(
            "--stage", type=Stage.parse, default=Stage.ALL,
            help="comma-separated list of processing stages (default: all)"
        )
        ap.add_argument(
            "--no-download", dest="download", action="store_false",
            help="do not download missing raw6 files from MAST",
        )
        ap.add_argument(
            "--recreate", action="store_true",
            help="erase and re-create any outputs that already exist",
        )
        ap.add_argument(
            "--verbose", "-v", action="count", default=0,
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
            "--aspect-dir", type=Path, default=None, metavar="DIR",
            help="find aspect files in this directory instead of"
            " <local_root>/aspect; does not need to be writable"
        )
        ap.add_argument(
            "--eclipse-dir", type=Path, default=None, metavar="DIR",
            help="find and store per-eclipse data (raw6 and photonlist files)"
            " in this directory instead of <local_root>/eclipse; must be writable"
        )
        ap.add_argument(
            "--remote-eclipse-dir", type=Path, default=None, metavar="DIR",
            help="alternate location to check for per-eclipse data that is"
            " not already present in <eclipse_dir>; does not need to be writable"
        )
        ap.add_argument(
            "--depth", type=parse_depth, default=DEFAULT_DEPTH,
            metavar="SECS",
            help="time (in seconds; fractions not allowed) to integrate"
            " over for each movie frame; use 'inf' or '-' for a single"
            " frame covering the entire observation leg"
        )
        ap.add_argument(
            "--aperture",
            type=parse_apertures, default=DEFAULT_APERTURES,
            metavar="ARCSEC[,ARCSEC,...]",
            help="comma-separated list of aperture sizes (in arcseconds)"
            " to use to compute photometry"
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
        ap.add_argument(
            "--print-tracebacks", action="store_true",
            help="print tracebacks for Python exceptions"
        )

        ns = ap.parse_args(args)

        if (aspect_dir := ns.aspect_dir) is None:
            aspect_dir = ns.catalog_root / "aspect"

        if (eclipse_dir := ns.eclipse_dir) is None:
            eclipse_dir = ns.catalog_root / "eclipse"

        return cls(
            eclipses = ns.eclipse,
            stages = ns.stage,
            bands = ns.band,
            depth = ns.depth,
            apertures = ns.aperture,
            download = ns.download,
            recreate = ns.recreate,
            verbose = ns.verbose,
            parallel = ns.parallel,
            chunk_size = ns.chunk_size,
            share_memory = ns.share_memory,
            catalog_root = ns.catalog_root,
            aspect_dir = aspect_dir,
            eclipse_dir = eclipse_dir,
            remote_eclipse_dir = ns.remote_eclipse_dir,
            print_tracebacks = ns.print_tracebacks,
        )


def report_uncaught_exception(
    e: Exception,
    eclipse: int,
    print_tracebacks: bool,
) -> None:
    if print_tracebacks:
        sys.stderr.write(
            f"\nException while processing eclipse {eclipse}:\n"
        )
        import traceback
        traceback.print_exc()
        sys.stderr.write("\n")
        return

    if isinstance(e, OSError):
        # Convert an OSError to a human-readable message the way I think
        # it should be done, which is different from the way the Python
        # devs think it should be done.  Also deals with OSError objects
        # created by e.g. `raise FileNotFoundError(path)`, which will be
        # structured completely differently than the equivalent thrown by
        # a failing `open(path)`.  (We care because pyarrow does this.)

        if e.strerror:
            # if e.strerror is set, we can safely assume this exception
            # was created by CPython core and e.filename/e.filename2
            # will also be set when relevant
            if e.filename and e.filename2:
                msg = f"{e.filename!r} -> {e.filename2!r}: {e.strerror}"
            elif e.filename:
                msg = f"{e.filename!r}: {e.strerror}"
            else:
                msg = e.strerror
        else:
            import errno
            from os import strerror

            # we assume this list to be comprehensive:
            # https://docs.python.org/3.10/library/exceptions.html#os-exceptions
            # ConnectionError is an intermediate base class and we'd do the
            # same thing for it that we do for the catch-all case at the
            # bottom, so it's excluded
            match e:
                # Unambiguous cases
                case ChildProcessError():
                    err = strerror(errno.ECHILD)
                case ConnectionAbortedError():
                    err = strerror(errno.ECONNABORTED)
                case ConnectionRefusedError():
                    err = strerror(errno.ECONNREFUSED)
                case ConnectionResetError():
                    err = strerror(errno.ECONNRESET)
                case FileExistsError():
                    err = strerror(errno.EEXIST)
                case FileNotFoundError():
                    err = strerror(errno.ENOENT)
                case InterruptedError():
                    err = strerror(errno.EINTR)
                case IsADirectoryError():
                    err = strerror(errno.EISDIR)
                case NotADirectoryError():
                    err = strerror(errno.ENOTDIR)
                case ProcessLookupError():
                    err = strerror(errno.ESRCH)
                case TimeoutError():
                    err = strerror(errno.ETIMEDOUT)

                # Ambiguous cases
                # BlockingIOError can mean either EWOULDBLOCK or EINPROGRESS.
                # (The Python docs also list their respective aliases EAGAIN
                # and EALREADY.)  Either condition reaching this function is
                # almost surely a bug, so we just pick one.
                case BlockingIOError():
                    err = strerror(errno.EWOULDBLOCK)

                # BrokenPipeError can mean either EPIPE or ESHUTDOWN.
                # The latter is more likely to come up in this application.
                case BrokenPipeError():
                    err = strerror(errno.ESHUTDOWN)

                # PermissionError also covers EPERM, but EPERM almost always
                # means "you tried to do something that only root can do"
                # (e.g. set the system's idea of "wall clock" time), which
                # should never come up in this application.
                case PermissionError():
                    err = strerror(errno.EACCES)

                # Subclasses that were added after Python 3.10 will wind up
                # here.  Hopefully the type's name will be helpful.
                case other:
                    err = type(other).__name__
                    if err == "OSError":
                        err = "system error, details lost"

            # guess that any string payload of the exception is a relevant
            # file name
            if e.args:
                msg = (
                    " -> ".join(repr(a) for a in e.args)
                    + ": "
                    + err
                )
            else:
                msg = err
    else:
        # Other exception types are much less troublesome, but we do have
        # to watch out for exceptions that stringify to ''.
        msg = str(e)
        if not msg:
            msg = type(e).__name__

    sys.stderr.write(f"glcat: eclipse {eclipse}: error: {msg}\n")


def process_eclipse(eclipse: int, options: CLIOptions):
    e_warn, e_error, postcsp = check_eclipse(eclipse, options.aspect_dir)
    for warning in e_warn:
        sys.stderr.write(f"glcat: eclipse {eclipse}: warning: {warning}\n")
    for err in e_error:
        sys.stderr.write(f"glcat: eclipse {eclipse}: error: {err}\n")
    if postcsp:
        print(f"Eclipse {eclipse} is post CSP.")
    if e_error:
        raise RuntimeError("skipped due to metadata issues\n")

    if options.stages & Stage.DOWNLOAD:
        stages.download_raw(
            eclipse,
            eclipse_dir=options.eclipse_dir,
            remote_eclipse_dir=options.remote_eclipse_dir,
            bands=options.bands,
            download=options.download,
            recreate=options.recreate,
            verbose=options.verbose,
        )
    if options.stages & Stage.BASE_PHOTOMETRY:
        stages.base_photometry(
            eclipse,
            eclipse_dir=options.eclipse_dir,
            remote_eclipse_dir=options.remote_eclipse_dir,
            aspect_dir=options.aspect_dir,
            bands=options.bands,
            depth=options.depth,
            aperture_sizes=options.apertures,
            recreate=options.recreate,
            verbose=options.verbose,
            parallel=options.parallel,
            chunk_size=options.chunk_size,
            share_memory=options.share_memory
        )
    if options.stages & Stage.FORCED_PHOTOMETRY:
        stages.forced_photometry(
            eclipse,
            eclipse_dir=options.eclipse_dir,
            remote_eclipse_dir=options.remote_eclipse_dir,
            aspect_dir=options.aspect_dir,
            bands=options.bands,
            depth=options.depth,
            aperture_sizes=options.apertures,
            recreate=options.recreate,
            verbose=options.verbose,
            parallel=options.parallel,
            chunk_size=options.chunk_size,
            share_memory=options.share_memory
        )
    if options.stages & Stage.BAND_CATALOG:
        stages.band_catalog(
            eclipse,
            eclipse_dir=options.eclipse_dir,
            aspect_dir=options.aspect_dir,
            bands=options.bands,
            recreate=options.recreate,
            verbose=options.verbose,
        )
    if options.stages & Stage.MERGED_CATALOG:
        stages.merged_catalog(
            eclipse,
            eclipse_dir=options.eclipse_dir,
            aspect_dir=options.aspect_dir,
            recreate=options.recreate,
            verbose=options.verbose,
        )


def main() -> int:
    options = CLIOptions.from_args()

    # TODO: cross-stage/cross-eclipse parallelism?
    # spin up a multiprocessing Pool and farm out invocations
    # of process_eclipse.  requires figuring out how to split
    # options.parallel into eclipse-level and chunk-level parallelism.
    status = 0
    for eclipse in options.eclipses:
        try:
            process_eclipse(eclipse, options)
        except Exception as e:
            status = 1
            report_uncaught_exception(e, eclipse, options.print_tracebacks)
    return status


if __name__ == "__main__":
    sys.exit(main())
