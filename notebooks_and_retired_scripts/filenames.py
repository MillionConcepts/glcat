""" utilities for making filenames for glcat """
import sys
sys.path.append('/home/bekah/gphoton_working/gPhoton')
from gPhoton.types import Pathlike, GalexBand
from Typing import Literal


def get_normal_frame_files(frame):
    """ edit to add final file names """
    frame_info = {}
    frame_info['1s_dose'] =
    frame_info['xy_list'] =
    frame_info['ra'] =
    frame_info['dec'] =

    return frame_info


def get_slew_frame_files(frame):
    """ edit to add final file names """
    frame_info = {}
    frame_info['1s_dose'] =
    frame_info['smoothed_1s_dose'] =
    frame_info['streaks'] =
    frame_info['stacked'] =
    frame_info['ra'] =
    frame_info['dec'] =
    frame_info['snippet'] =
    #frame_info['h'] =
    #frame_info['w'] = # ra and dec guess for h and w

    return frame_info


def get_main_files(
    eclipse: int,
    band: GalexBand = "NUV",
    compression: Literal["none", "gzip", "rice"] = "gzip",
    root: Pathlike = "data",
    **kwargs,
) -> dict[str, str]:
    """
    modified from gphoton eclipse_to_paths.
    """
    root = "data" if root is None else root

    zpad = str(eclipse).zfill(5)

    eclipse_path = f"{root}/e{zpad}/"

    eclipse_base = f"{eclipse_path}e{zpad}"

    ext = {"gzip": ".fits.gz", "none": ".fits", "rice": ".fits"}[compression]
    comp = {"gzip": "g", "none": "u", "rice": "r"}[compression]
    # automatically select direct, leaving in in case that changes
    mode = {"direct": "d", "grism": "g", "opaque": "o"}['direct']

    prefix = f"{eclipse_base}-{band[0].lower()}{mode}"
    file_dict = {
        "raw6": f"{prefix}-raw6.fits.gz",
        "photonfile": f"{prefix}.parquet",
        "old_aspect_parq":

    }

    return file_dict



#********************************************************************
# functions may note use (went normal vs slew route)

def get_backplane_filename(eclipse, t, expt, root):
    """ t is timestamp in aspect table associated with backplane"""

    depth = None if depth is None else f"t{str(depth).zfill(4)}"

    return backplane_file


def get_xylist_filename(eclipse, t, expt, root):

    depth = None if depth is None else f"t{str(depth).zfill(4)}"

    return xylist_file


def get_wcs_filename(eclipse, t, expt, root):
    depth = None if depth is None else f"t{str(depth).zfill(4)}"

    return wcs_file