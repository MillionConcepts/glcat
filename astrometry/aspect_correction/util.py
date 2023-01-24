# helper functions for aspect refinement

from pyarrow import parquet
from astropy.io import fits
import numpy as np
import os
import pandas as pd
import subprocess
import io


def make_refined_aspect_table(file_names):
    """makes pandas df with wcs info for all frames, including failed ones."""
    # making an empty dataframe with column names, could pull col names from first
    # returned df but would get messy if first few frames failed
    asp_df = pd.DataFrame(columns=['frame', 'crpix0', 'crpix1', 'crval0', 'ra_tangent',
                                   'dec_tangent', 'pixx_tangent', 'pixy_tangent',
                                   'imagew', 'imageh', 'cd11', 'cd12', 'cd21', 'cd22',
                                   'det', 'parity', 'pixscale', 'orientation',
                                   'ra_center', 'dec_center', 'orientation_center',
                                   'ra_center_h', 'ra_center_m', 'ra_center_s',
                                   'dec_center_sign', 'dec_center_d', 'dec_center_m',
                                   'dec_center_s', 'ra_center_hms', 'dec_center_dms',
                                   'ra_center_merc', 'dec_center_merc', 'fieldarea',
                                   'fieldw', 'fieldh', 'fieldunits', 'decmin', 'decmax',
                                   'ramin', 'ramax', 'ra_min_merc', 'ra_max_merc',
                                   'dec_min_merc', 'dec_max_merc', 'merc_diff',
                                   'merczoom', 'failed_flag'])
    for frame, wcs_path in file_names["output_wcs"]:
        if os.path.exists(wcs_path):
            frame_wcs = get_aspect_from_wcsinfo(wcs_path)
            frame_wcs["frame"] = frame
            asp_df = pd.concat([asp_df, frame_wcs])
        if not os.path.exists(wcs_path):
            print("This frame failed via xylist, marking as failed in aspect df.")
            frame_no_wcs = pd.DataFrame(columns=['frame', 'failed_flag'])
            frame_no_wcs["frame"] = frame
            frame_no_wcs["failed_flag"] = True
            asp_df = pd.concat([asp_df, frame_no_wcs])

    print("The list of failed frames is:")
    print(asp_df[asp_df['failed_flag']==True]['frame'])

    return asp_df.set_index('frame')


def make_file_names(opt):
    """ make dict of file names for output files produced in aspect refinement
    (can also be used to delete files later, because we don't need to keep most
     of these). not all of these names will be used for every run time, but
      strings don't take up a lot of room. """

    #TODO: don't hardcode all the path names :)

    main_path = opt['output_dir']
    aspect_solns = f'{main_path}aspect_solns/'
    xylist_folder = f'{main_path}xylists_lower_thresh/e09869/'
    astrometry_temp = f'{main_path}astrometry_temp/'

    eclipse = str(opt['eclipse']).zfill(5)

    # main refined aspect solution
    if opt['threshold'] and opt['star_size'] is None:
        asp_df = f"{aspect_solns}e{eclipse}_aspect_soln_{opt['expt']}s"
    else:
        asp_df = f"{aspect_solns}e{eclipse}_aspect_soln_{opt['expt']}s_thresh{opt['threshold']}" \
                 f"_size{opt['star_size']}_dose_gaia_tycho"

    fits_table = []
    frame_path = []
    output_wcs = []
    ver_wcs = []
    xylist = []

    if not opt['dose']:
        image_path = f'/home/bekah/glcat/astrometry/e{eclipse}/'
        for i in range(opt['num_frames']):
            # for feeding xylist to astrometry.net
            # think there's a fancier list comprehension way to do this but can't google bc on a plane
            fits_table.append(f'/home/bekah/glcat/astrometry/aspect_correction/frame{i}.xyls')
            # paths to movie frame images, largely used for an image or verification run
            frame_padded = str(i).zfill(4)
            frame_path.append(str(image_path + f"e{eclipse}-nd-{opt['expt']}s-"
                                               f"f{frame_padded}-rice.fits")) # 0-
            # for each frame, resulting new wcs
            output_wcs.append((i, f"/home/bekah/glcat/astrometry/aspect_correction/frame{i}.wcs")) # tuple of frame and wcs path
            # verified wcs list
            ver_wcs.append(f"/home/bekah/glcat/astrometry/aspect_correction/e21442/e21442-"
                           f"nd-1s-f{frame_padded}-rice.wcs")
        # may not be used, for output of a verification run
        verified_df = f"{eclipse}_aspect_soln_{opt['expt']}s_thresh{opt['threshold']}_" \
                      f"size{opt['star_size']}_verified"
        # main refined aspect solution
        if opt['threshold'] and opt['star_size'] is None:
            asp_df = f"{eclipse}_aspect_soln_{opt['expt']}s"
        else:
            asp_df = f"{eclipse}_aspect_soln_{opt['expt']}s_thresh{opt['threshold']}_" \
                     f"size{opt['star_size']}_noMinus1"

    if opt['dose']:
        image_path = f'/home/bekah/glcat/astrometry/e{eclipse}/dose_t2_selection_2/'
        for i in range(opt['num_frames']):
            # for feeding xylist to astrometry.net
            # think there's a fancier list comprehension way to do this but can't google bc on a plane
            fits_table.append(f'{xylist_folder}frame{i}.xyls')
            # paths to movie frame images, largely used for an image or verification run
            frame_padded = str(i).zfill(5)
            time_padded = str(opt['expt']).zfill(4)
            eclipse = str(eclipse).zfill(5)
            frame_path.append(str(image_path + f"e{eclipse}-nd-t{time_padded}-b00-"
                                               f"f{frame_padded}-g_dose.fits")) # 0-
            # for each frame, resulting new wcs
            output_wcs.append(f"{astrometry_temp}frame{i}.wcs")
            # verified wcs list
            ver_wcs.append(f"{astrometry_temp}e21442-nd-1s-f{frame_padded}-rice.wcs")

        # may not be used, for output of a verification run
        verified_df = f"{astrometry_temp}e{eclipse}_aspect_soln_{opt['expt']}s_" \
                      f"thresh{opt['threshold']}_size{opt['star_size']}_verified"

    file_names = {"asp_df": asp_df, "verified_df": verified_df, "xylist": fits_table,
                  "frame_path": frame_path, "output_wcs": output_wcs, "ver_wcs": ver_wcs,
                  "astrometry_temp": astrometry_temp}
    return file_names


def get_ra_dec(eclipse):
    """ Loading aspect table, set the ra and dec to 0 in the pipeline to get the movie,
     so we have to work around that here to get an initial ra / dec guess for
      astrometry.net. """
    eclipse = 9869
    parq = parquet.read_table('/home/bekah/gphoton_working/gPhoton/aspect/aspect.parquet')
    aspect = parq.to_pandas()
    ra = np.mean(aspect[aspect["eclipse"] == eclipse]["ra"])
    dec = np.mean(aspect[aspect["eclipse"] == eclipse]["dec"])
    return ra, dec


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


def zero_flag_and_edge(cnt, flag, edge):
    """ set flag and edge pixels to zero within count array, used to prevent
    fake sources. """
    cnt[~np.isfinite(cnt)] = 0
    cnt[np.nonzero(flag)] = 0
    cnt[np.nonzero(edge)] = 0
    return cnt


#---- OLD FUNCTIONS -----

def clean_up(file_names, i):
    #cmd = f"rm glcat_tests/astrometry_temp/*"
    #subprocess.call(cmd, shell=True)
    return None

def get_aspect_from_wcs(wcs_path):
    """To get the new ra and dec for each frame, open the output wcs file."""
    #TODO: it would be nice to not have to close and open so many files, also check
    # that these files actually each close after use
    frame_wcs = fits.open(wcs_path)
    new_center_ra = frame_wcs[0].header["CRVAL1"]
    new_center_dec = frame_wcs[0].header["CRVAL2"]
    return new_center_ra, new_center_dec

