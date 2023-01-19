""" Main pipeline for aspect refinement on individual frames (short time period, not projected)
    of eclipses. Sources are detected using DAO star finder, source list + intiial ra / dec
    are fed to astrometry.net when using an "xylist" run.

    Uses astrometry.net (https://astrometry.net/) takes astronomical images and uses "catalogs"
     to solve for refined aspect solutions. """

from pyarrow import parquet
from astropy.table import Table
import sys
import numpy as np
import fitsio
import pandas as pd
from astropy.io import fits
from astropy import wcs
import numpy.ma as ma
import subprocess
import os.path
from typing import Optional
from xylist_runs import run_xylist
from image_runs import run_image_frame, run_image, run_verification
from util import make_file_names


def execute_refiner(
                    eclipse,
                    expt,
                    num_frames,
                    run_type,
                    dose: Optional[bool] = False,
                    threshold: Optional[float] = None,
                    star_size: Optional[int] = None,
                    verify: Optional[bool] = False,
                    output_directory: Optional[str] = "/home/bekah/glcat_tests"
                    ):
    """
    Args:
        eclipse: for calling files
        expt: number of seconds per frame used in gphoton2, for file name
        num_frames: can come up with a better way to do this, currently just
         for knowing how many times to run loop2
        run_type: can be "image" or "xylist". xylist runs DAOstarfinder
         and feeds list of sources to astrometry.net, image feeds the image
         directly to astrometry.net. image tends to take longer.
        dose: are we using dose maps?
        threshold: threshold for star detection used by dao, optional
        star_size: what size of star dao is looking for, will probably
         vary with the frame length used in the run
        verify: runs a second run of astrometry.net for "verification" that
         uses outputs of previous run. may take a while, but hopefully improves
         output.
        output_directory: where to save the files
     """
    print(f"Executing refiner for eclipse {eclipse} using astrometry.net")
    print(f"there are {num_frames} frames for the eclipse.")

    # setting options and padding eclipse num
    opt = {
        'threshold': threshold,
        'star_size': star_size,
        'dose': dose,
        'verify': verify,
        'output_dir': output_directory
    }

    print(f"Making file name dictionary. File names will be saved to: {output_directory}")
    file_names = make_file_names(str(eclipse).zfill(5), expt, num_frames, run_type, threshold, star_size)

    # if you want to rerun the aspect correction to hopefully "refine" the solution, ends run
    # after verification
    if verify:
        print("Running verification on previous aspect refinement.")
        verified_df = run_verification(eclipse, expt, num_frames, file_names, **opt)
        verified_df.to_csv(file_names["verified_df"])
        return

    if run_type == "image":
        # useful for testing the astrometry.net "in house" source extraction's performance
        print("Running astrometry.net directly on a frame images.")
        asp_df = run_image(eclipse, expt, num_frames, file_names)
    elif run_type == "xylist":
        print("Running astrometry.net on xylists derived from frame images.")
        asp_df = run_xylist(eclipse, expt, num_frames, file_names, dose, **opt)
    print("writing aspect solutions to csv.")
    asp_df.to_csv(file_names["asp_df"])

    return




