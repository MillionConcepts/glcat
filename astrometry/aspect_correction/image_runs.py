# Functions for running astrometry.net using images and no DAO star finder

from astropy.io import fits
import subprocess
import os.path
from util import get_ra_dec, make_refined_aspect_table


def run_image(eclipse, expt, num_frames, file_names):
    ra, dec = get_ra_dec(eclipse)
    for i in range(num_frames):
        frame_path = file_names["frame_path"][i]
        print("running astrometry.net on IMAGE from individual frame of movie")
        print(frame_path)
        movie = fits.open(frame_path)
        image_width = movie[1].header['NAXIS1']
        image_height = movie[1].header['NAXIS2']
        crpix_x = movie[1].header['CRPIX1']   # RA -- TAN
        crpix_y = movie[1].header['CRPIX2']   # DEC -- TAN
        cmd = f"solve-field --overwrite --sigma 0.3 -y -3 {ra} -4 {dec} " \
              f"-w {image_width} -e {image_height} -L 1.2 -H 1.8 -u app " \
              f"  --crpix-x {crpix_x} --crpix-y {crpix_y} --radius 5" \
              f" --cpulimit 500 --extension 1 {frame_path} "
        # --sigma 0.5,  --no-tweak
        subprocess.call(cmd, shell=True)
    asp_df = make_refined_aspect_table(file_names)
    return asp_df


def run_image_frame(eclipse, expt, i, file_names):
    #  for frame that fails w xy-list
    ra, dec = get_ra_dec(eclipse)
    frame_path = file_names["frame_path"][i]
    print("running astrometry.net on individual frame of movie")
    print(frame_path)
    movie = fits.open(frame_path)
    image_width = movie[0].header['NAXIS1'] # used to be index 1 for all of these when not dose maps
    image_height = movie[0].header['NAXIS2']
    crpix_x = movie[0].header['CRPIX1']   # RA -- TAN
    crpix_y = movie[0].header['CRPIX2']   # DEC -- TAN
    cmd = f"solve-field --overwrite --sigma 0.2 -y -3 {ra} -4 {dec} " \
              f"-w {image_width} -e {image_height} --scale-units arcsecperpix" \
              f" --scale-low 1.4 --scale-high 1.6" \
              f"  --crpix-x {crpix_x} --crpix-y {crpix_y} --radius 5" \
              f"  --cpulimit 600 --no-plots --extension 1 {frame_path} "
# --sigma 0.5,  --no-tweak
    subprocess.call(cmd, shell=True)
    wcs_path = file_names["output_wcs"][i]
    return wcs_path


def run_verification(eclipse, expt, num_frames, file_names, **opt):
    """Use the verification mode of astrometry.net to check an existing wcs file."""
    aspect_solution = {}
    ra, dec = get_ra_dec(eclipse)
    for frame in range(num_frames):
        frame_path = file_names["frame_path"][frame]
        print("Running astrometry.net on image of movie frame.")
        movie = fits.open(frame_path)
        image_width = movie[1].header['NAXIS1']
        image_height = movie[1].header['NAXIS2']
        crpix_x = movie[1].header['CRPIX1']  # RA -- TAN
        crpix_y = movie[1].header['CRPIX2']  # DEC -- TAN
        old_wcs = file_names["output_wcs"][frame]
        if os.path.exists(old_wcs):
            cmd = f"solve-field -v --overwrite --verify {old_wcs} -y -3 {ra} -4 {dec} " \
                  f"-w {image_width} -e {image_height} -L 1.2 -H 1.8 -u app " \
                  f"  --crpix-x {crpix_x} --crpix-y {crpix_y} --radius 5" \
                  f" --cpulimit 500 --extension 1 {frame_path} "
            # --sigma 0.5,  --no-tweak
            subprocess.call(cmd, shell=True)
    asp_df = make_refined_aspect_table(file_names)
    return asp_df