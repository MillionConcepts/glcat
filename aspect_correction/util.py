# helper functions for aspect refinement

from pyarrow import parquet
import numpy as np
import os
import pandas as pd
import subprocess
import io
from astropy.io import fits
from typing import Optional


def crop_xylist(xylist, xylist_cropped, x_dim, y_dim):

    if x_dim != y_dim:
        print("XYlist cannot be cropped because image is not a square.")
        return
    elif x_dim == y_dim:
        xy_list = fits.open(xylist)
        tbdata = xy_list[1].data
        # cropping to between 700 and 2500 in the x and y directions
        selected = (tbdata['X'] > 200) & (tbdata['X'] < 3000) & (tbdata['Y'] > 200) & (tbdata['Y'] < 3000)
        newtbdata = tbdata[selected]
        # resetting so that the origin of cropped area is 0
        hdu = fits.BinTableHDU(data=newtbdata)
        hdu.data['X'] = hdu.data['X'] - 200
        hdu.data['Y'] = hdu.data['Y'] - 200
        hdu.writeto(xylist_cropped, overwrite=True)
        # 700 subtracted twice from both axis, so it's 1400
        new_dims = (3200-400, 3200-400)

    return new_dims


def make_refined_aspect_table(file_names, crop: Optional[bool] = False):
    """makes pandas df with wcs info for all frames, including failed ones."""
    # making an empty dataframe with column names, could pull col names from first
    # returned df but would get messy if first few frames failed
    asp_df = pd.DataFrame(columns=['frame', 'crpix0', 'crpix1', 'crval1', 'crval0', 'ra_tangent',
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
    if crop:
        wcs_list = file_names["output_wcs_cropped"]
    else:
        wcs_list = file_names["output_wcs"]

    for frame, wcs_path in wcs_list:
        if os.path.exists(wcs_path):
            frame_wcs = get_aspect_from_wcsinfo(wcs_path)
            frame_wcs["frame"] = frame
            asp_df = pd.concat([asp_df, frame_wcs])
        if not os.path.exists(wcs_path):
            frame_no_wcs = pd.DataFrame(columns=['frame', 'failed_flag'])
            frame_no_wcs["frame"] = frame
            frame_no_wcs["failed_flag"] = True
            asp_df = pd.concat([asp_df, frame_no_wcs])
    asp_df = asp_df.set_index('frame')
    asp_df = asp_df.astype(
        {'crpix0': 'int64',
         'crpix1': 'int64',
         'crval0': 'float64',
         'ra_tangent': 'float64',
         'dec_tangent': 'float64',
         'pixx_tangent': 'int64',
         'pixy_tangent': 'int64',
         'imagew': 'int64',
         'imageh': 'int64',
         'cd11': 'float64',
         'cd12': 'float64',
         'cd21': 'float64',
         'cd22': 'float64',
         'det': 'float64',
         'parity': 'int64',
         'pixscale': 'float64',
         'orientation': 'float64',
         'ra_center': 'float64',
         'dec_center': 'float64',
         'orientation_center': 'float64',
         'ra_center_h': 'int64',
         'ra_center_m': 'int64',
         'ra_center_s': 'float64',
         'dec_center_sign': 'int64',
         'dec_center_d': 'int64',
         'dec_center_m': 'int64',
         'dec_center_s': 'float64',
         'ra_center_hms': 'O',
         'dec_center_dms': 'O',
         'ra_center_merc': 'float64',
         'dec_center_merc': 'float64',
         'fieldarea': 'float64',
         'fieldw': 'float64',
         'fieldh': 'float64',
         'fieldunits': 'O',
         'decmin': 'float64',
         'decmax': 'float64',
         'ramin': 'float64',
         'ramax': 'float64',
         'ra_min_merc': 'float64',
         'ra_max_merc': 'float64',
         'dec_min_merc': 'float64',
         'dec_max_merc': 'float64',
         'merc_diff': 'float64',
         'merczoom': 'int64',
         'failed_flag': 'float64',
         'crval1': 'float64'})
    # adds rows for failed frames
    #asp_df = asp_df.reindex(range(1641), method=None)
    #print("The list of failed frames is:")
    #null_list = asp_df.index[asp_df.isnull().all(1)].to_list()
    #null_df = {'failed_frames': null_list}
    #null_df = pd.DataFrame(null_df)
    #print(null_list)
    #null_df.to_csv(file_names["null_list"])
    # linearly interpolates
    asp_df = asp_df.interpolate(method='linear', limit_direction='both', axis=0)
    return asp_df


def make_file_names(opt):
    """ make dict of file names for output files produced in aspect refinement
    (can also be used to delete files later, because we don't need to keep most
    of these). not all of these names will be used for every run time, but
    strings don't take up a lot of room.
    assumes you have folders in your main director called 'aspect_solns' and
    'xylists' and 'astrometry_temp' """

    # directories
    main_path = opt['output_dir']
    aspect_solns = f'{main_path}aspect_solns/'
     #f"{main_path}xylists_9871/e09871/" #f'{main_path}xylists_lower_thresh/e09869_0.8_3/'
    astrometry_temp = f'{main_path}astrometry_temp/'
    eclipse = str(opt['eclipse']).zfill(5)
    xylist_folder = f"/home/bekah/gphoton_working/test_data/e{eclipse}/"
    dose_image_path = f'/home/bekah/glcat/astrometry/e{eclipse}/dose_t2_selection_2/'
    image_path = f'/home/bekah/glcat/astrometry/e{eclipse}/'
    # main refined aspect solution
    #if opt['threshold'] and opt['star_size'] is None:
    #    if opt["crop"]:
    #      asp_df = f"{aspect_solns}e{eclipse}_aspect_soln_{opt['expt']}_s_cropped"
    #     else:
    #        asp_df = f"{aspect_solns}e{eclipse}_aspect_soln_{opt['expt']}_s"
    #else:
    #    if opt["crop"]:
    #        asp_df = f"{aspect_solns}e{eclipse}_aspect_soln_{opt['expt']}s_thresh" \
    #             f"{opt['threshold']}_size{opt['star_size']}_cropped"
    #    else:
    #        asp_df = f"{aspect_solns}e{eclipse}_aspect_soln_{opt['expt']}s_thresh" \
    #         f"{opt['threshold']}_size{opt['star_size']}_dose_gaia_tycho"
    # image paths
    os.mkdir(f"/home/bekah/gphoton_working/astrom_test_data/e{eclipse}/")
    asp_df = f"/home/bekah/gphoton_working/astrom_test_data/e{eclipse}/astrometry_xy_{opt['expt']}s_{eclipse}"
    frame_path = []
    if not opt['dose']:
        for i in range(opt['num_frames']):
            # paths to movie frame images, largely used for an image or verification run
            frame_padded = str(i).zfill(4)
            frame_path.append(str(image_path + f"e{eclipse}-nd-{opt['expt']}s-"
                                               f"f{frame_padded}-rice.fits")) # 0-
    if opt['dose']:
        for i in range(opt['num_frames']):
            # paths to movie frame images, largely used for an image or verification run
            frame_padded = str(i).zfill(5)
            time_padded = str(opt['expt']).zfill(4)
            frame_path.append(str(dose_image_path + f"e{eclipse}-nd-t{time_padded}-b00-"
                                               f"f{frame_padded}-g_dose.fits")) # 0-
    # wcs and xylist paths
    output_wcs = []
    output_wcs_cropped = []
    ver_wcs = []
    fits_table = []
    fits_table_cropped = []
    for i in range(opt['num_frames']):
        frame_padded = str(i).zfill(5)
        # FRAME EXTENSION (CHANGE)
        frame_extension ="" #f"_{opt['expt']}s_e{eclipse}"
        fits_table.append(f'{xylist_folder}frame{i}{frame_extension}.xyls')
        # for if we want to crop an area out of xylist by coordinates
        fits_table_cropped.append(f'{xylist_folder}frame{i}_crop.xyls')
        # for each frame, resulting new wcs
        output_wcs.append((i, f"{astrometry_temp}frame{i}{frame_extension}.wcs"))
        # for cropped xylists wcs
        output_wcs.append((i, f"{astrometry_temp}frame{i}{frame_extension}.wcs"))
        # verified wcs list
        ver_wcs.append(f"{astrometry_temp}e21442-nd-1s-f{frame_padded}-rice.wcs")
    # may not be used, for output of a verification run
    verified_df = f"{astrometry_temp}e{eclipse}_aspect_soln_{opt['expt']}s_" \
                      f"thresh{opt['threshold']}_size{opt['star_size']}_verified"
    null_list = f"{aspect_solns}e{eclipse}_null_frames_{opt['expt']}s.csv"
    file_names = {"asp_df": asp_df, "verified_df": verified_df, "xylist": fits_table,
                  "xylist_cropped": fits_table_cropped,
                  "frame_path": frame_path, "output_wcs": output_wcs,
                  "output_wcs_cropped": output_wcs_cropped, "ver_wcs": ver_wcs,
                  "astrometry_temp": astrometry_temp, "null_list": null_list}
    return file_names


def get_ra_dec(eclipse):
    """ Loading aspect table, set the ra and dec to 0 in the pipeline to get the movie,
     so we have to work around that here to get an initial ra / dec guess for
      astrometry.net. """
    # have to change parquet table location for another user
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



