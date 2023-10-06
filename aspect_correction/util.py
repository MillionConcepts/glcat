# helper functions for aspect refinement

from pyarrow import parquet
import numpy as np
import os
import pandas as pd
import subprocess
import io
from astropy.io import fits
from typing import Optional


def get_aspect_from_wcsinfo(wcs_path):
    """returns 1 row pd df with columns that are wcs values read from the
    astrometry.net output wcs file read w/ astrometry.net program wcsinfo,
     we use wcsinfo instead of reading the wcs header because it returns more
     information using other functions within astrometry.net (ie could
     recreate wcsinfo ourselves but it's a lot of work and written in c)."""
    output = subprocess.check_output(f"wcsinfo {wcs_path}",
                                     stderr=subprocess.STDOUT,
                                     shell=True)
    output = output.decode()  # bc was bytes not str
    output = "name value\n"+output
    frame_wcs = pd.read_csv(io.StringIO(output), sep=' ')
    frame_wcs = frame_wcs.set_index('name')
    frame_wcs = frame_wcs.transpose()
    return frame_wcs


def get_quality_metrics_from_wcs(wcs_path):
    """ read quality metrics from wcs fits header """
    hdu = fits.open(wcs_path)
    hdu[0].header
    return


def zero_flag_and_edge(cnt, flag, edge):
    """ set flag and edge pixels to zero within count array, used to prevent
    fake sources. """
    cnt[~np.isfinite(cnt)] = 0
    cnt[np.nonzero(flag)] = 0
    cnt[np.nonzero(edge)] = 0
    return cnt


# *********************************************************************************
# below function (execute combine) requires functions which have been retired
# so it will have to be updated to be used
# from util import make_refined_aspect_table, make_file_names
# from typing import Optional

def execute_combine(
                    eclipse,
                    expt,
                    num_frames,
                    run_type,
                    dose: Optional[bool] = False,
                    threshold: Optional[float] = None,
                    star_size: Optional[int] = None,
                    output_directory: Optional[str] = "/home/bekah/glcat_tests/"
                    ):
    opt = {
        'eclipse': eclipse,
        'expt': expt,
        'num_frames': num_frames,
        'threshold': threshold,
        'star_size': star_size,
        'dose': dose,
        'run_type': run_type,
        'output_dir': output_directory
    }
    #file_names = make_file_names(opt)
    #asp_table = make_refined_aspect_table(file_names)
    #asp_table.to_csv(file_names["asp_df"])
    return

