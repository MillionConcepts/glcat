#  Uses movies produced by gphoton2 to solve for new aspect solution using astrometry.net

# astrometry.net (https://astrometry.net/) takes astronomical images and uses "catalogs" to solve for
# refined aspect solutions

# with a known WCS header, solving happens a lot faster (not blindly guessing)
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


#TODO:
# 1)need to make it stop failing on random frames (why does this happen?), give longer lead time
# 2) use the 'dose' maps instead of projected stars

def execute_refiner(
                    eclipse,
                    expt,
                    num_frames,
                    run_type,
                    threshold: Optional[float] = None,
                    star_size: Optional[int] = None,
                    verify: Optional[bool] = False
                    ):
    """
    Main pipeline for aspect refinement on individual frames (short time period, not projected)
    of eclipses.
    Sources are detected using DAO star finder, source list + intiial ra / dec are fed to
    astrometry.net when using an "xylist" run.
     eclipse = for calling files
     expt = number of seconds per frame used in gphoton2, for file name
     num_frames = can come up with a better way to do this,
     currently just for knowing how many times to run loop2
     run_type = can be "image" or "xylist". xylist runs DAOstarfinder and feeds
     list of sources to astrometry.net, image feeds the image directly to
     astrometry.net. image tends to take longer.
     verify = runs a second run of astrometry.net for "verification" that uses outputs
     of previous run. may take a while, but hopefully improves output.
     """
    print(f"there are {num_frames} frames in the movie.")
    opt = {
        'threshold': threshold,
        'star_size': star_size,
    }
    if verify:
        print("running verification.")
        verified_df = run_verification(eclipse, expt, num_frames, **opt)
        verified_df.to_csv(f"{eclipse}_aspect_soln_{expt}s_thresh{threshold}_size{star_size}_verified")

    if run_type == "image":
        asp_df = run_image(eclipse, expt, num_frames)
    elif run_type == "xylist":
        asp_df = run_xylist(eclipse, expt, num_frames, **opt)
    print("writing aspect solutions to csv.")
    if threshold and star_size is None:
        asp_df.to_csv(f"{eclipse}_aspect_soln_{expt}s")
    else:
        asp_df.to_csv(f"{eclipse}_aspect_soln_{expt}s_thresh{threshold}_size{star_size}_wImageForFails")
    if verify:
        print("running verification.")
        verified_df = run_verification(eclipse, expt, num_frames, **opt)
        verified_df.to_csv(f"{eclipse}_aspect_soln_{expt}s_thresh{threshold}_size{star_size}_verified")
    return


def run_verification(eclipse, expt, num_frames, **opt):
    path = f'/home/bekah/glcat/astrometry/e{eclipse}/'

    aspect_solution = {}
    failed_frames = []

    path = f'/home/bekah/glcat/astrometry/e{eclipse}/'
    aspect_solution = {}
    ra, dec = get_ra_dec(eclipse)

    # for the first few frames that will inevitably fail
    new_center_ra = ra
    new_center_dec = dec

    for i in range(num_frames):

        frame_padded = str(i).zfill(4)
        frame_path = str(path + f"e{eclipse}-nd-{expt}s-f{frame_padded}-rice.fits")
        print("running astrometry.net on IMAGE from individual frame of movie")
        print(frame_path)

        movie = fits.open(frame_path)

        image_width = movie[1].header['NAXIS1']
        image_height = movie[1].header['NAXIS2']
        crpix_x = movie[1].header['CRPIX1']  # RA -- TAN
        crpix_y = movie[1].header['CRPIX2']  # DEC -- TAN

        old_wcs = f'/home/bekah/glcat/astrometry/aspect_correction/frame{i}.wcs'
        if os.path.exists(old_wcs):
            cmd = f"solve-field -v --overwrite --verify {old_wcs} -y -3 {ra} -4 {dec} " \
                  f"-w {image_width} -e {image_height} -L 1.2 -H 1.8 -u app " \
                  f"  --crpix-x {crpix_x} --crpix-y {crpix_y} --radius 5" \
                  f" --cpulimit 500 --extension 1 {frame_path} "
            # --sigma 0.5,  --no-tweak
            subprocess.call(cmd, shell=True)

        wcs_path = f"/home/bekah/glcat/astrometry/aspect_correction/e21442/e21442-nd-1s-f{frame_padded}-rice.wcs"

        if os.path.exists(wcs_path):
            new_center_ra, new_center_dec = get_aspect_from_wcs(wcs_path)
            aspect_solution[i] = [new_center_ra, new_center_dec]
    #     if not os.path.exists(wcs_path):2
    #         # try an alternative solution without DAO-produced XY lists
    #         run_image_frame(eclipse, expt, i)
    #         if os.path.exists(wcs_path):
    #             new_center_ra, new_center_dec = get_aspect_from_wcs(wcs_path)
    #             aspect_solution[i] = [new_center_ra, new_center_dec]
    #             # but use the previous ra and dec if that fails too
    #         if not os.path.exists(wcs_path):
    #             aspect_solution[i] = [new_center_ra, new_center_dec]
    #             print("This frame failed both via xylist and image.")
    #             failed_frames.append(i)
    #
    # print("The list of failed frames is:")
    # print(failed_frames)

    asp_df = pd.DataFrame(aspect_solution)

    return asp_df

def run_xylist(eclipse, expt, num_frames, **opt):

    path = f'/home/bekah/glcat/astrometry/e{eclipse}/'

    aspect_solution = {}
    failed_frames = []

    for i in range(num_frames):
        print(f"opening frame {i} fits file")
        frame_padded = str(i).zfill(4)
        frame_path = path+f"e{eclipse}-nd-{expt}s-f{frame_padded}-rice.fits" # -0, leg
        movie = fits.open(frame_path)
        print(f"masking flag and edge in frame {i}")
        masked_cnt = zero_flag_and_edge(movie[1].data[0], movie[2].data[0], movie[3].data[0])
        print(f"extracting sources from frame {i}")
        table_name = get_stars(masked_cnt, i, **opt)

        ra, dec = get_ra_dec(eclipse)
        if table_name is None and i == 0:
            print(f"No stars found in this frame {i}.")
            new_center_ra = ra
            new_center_dec = dec
            aspect_solution[i] = [ra, dec]
        if table_name is None:
            print(f"No stars found in this frame {i}.")
            aspect_solution[i] = [new_center_ra, new_center_dec]
        else:
            image_width = movie[1].header['NAXIS1']
            image_height = movie[1].header['NAXIS2']
            crpix_x = movie[1].header['CRPIX1']-1 # RA -- TAN
            crpix_y = movie[1].header['CRPIX2']-1 # DEC -- TAN
            print(f"crpix_x = {crpix_x}, crpix_y = {crpix_y}")
            print(f"image width = {image_width}, image height = {image_height}")

            run_astrometry_net(eclipse, expt, table_name, image_width, image_height, ra, dec, crpix_x, crpix_y)

            wcs_path = f"/home/bekah/glcat/astrometry/aspect_correction/frame{i}.wcs"
            if os.path.exists(wcs_path):
                new_center_ra, new_center_dec = get_aspect_from_wcs(wcs_path)
                aspect_solution[i] = [new_center_ra, new_center_dec]
            if not os.path.exists(wcs_path):
                # try an alternative solution without DAO-produced XY lists
                run_image_frame(eclipse, expt, i)
                if os.path.exists(wcs_path):
                    new_center_ra, new_center_dec = get_aspect_from_wcs(wcs_path)
                    aspect_solution[i] = [new_center_ra, new_center_dec]
                # but use the previous ra and dec if that fails too
                if not os.path.exists(wcs_path):
                    aspect_solution[i] = [new_center_ra, new_center_dec]
                    print("This frame failed both via xylist and image.")
                    failed_frames.append(i)

    print("The list of failed frames is:")
    print(failed_frames)

    asp_df = pd.DataFrame(aspect_solution)

    return asp_df


def run_image(eclipse, expt, num_frames):

    # TODO: need way to filter out images with only 0s
    #  (typically at beginnning of frames) bc it slows everything down

    path = f'/home/bekah/glcat/astrometry/e{eclipse}/'
    aspect_solution = {}
    ra, dec = get_ra_dec(eclipse)

    for i in range(num_frames):

        frame_padded = str(i).zfill(4)
        frame_path = str(path+f"e{eclipse}-nd-{expt}s-f{frame_padded}-rice.fits")
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

        wcs_path = f"/home/bekah/glcat/astrometry/aspect_correction/frame{i}.wcs"
        if os.path.exists(wcs_path):
            new_center_ra, new_center_dec = get_aspect_from_wcs(wcs_path)
            aspect_solution[i] = [new_center_ra, new_center_dec]

    asp_df = pd.DataFrame(aspect_solution)
    return asp_df


def run_image_frame(eclipse, expt, i):
    #  for frame that fails w xy-list

    path = f'/home/bekah/glcat/astrometry/e{eclipse}/'
    aspect_solution = {}
    ra, dec = get_ra_dec(eclipse)


    frame_padded = str(i).zfill(4)
    frame_path = str(path+f"e{eclipse}-nd-{expt}s-f{frame_padded}-rice.fits")
    print("running astrometry.net on individual frame of movie")
    print(frame_path)

    movie = fits.open(frame_path)

    image_width = movie[1].header['NAXIS1']
    image_height = movie[1].header['NAXIS2']
    crpix_x = movie[1].header['CRPIX1']   # RA -- TAN
    crpix_y = movie[1].header['CRPIX2']   # DEC -- TAN

    cmd = f"solve-field --overwrite --sigma 0.2 -y -3 {ra} -4 {dec} " \
              f"-w {image_width} -e {image_height} -L 1 -H 1.8 -u app " \
              f"  --crpix-x {crpix_x} --crpix-y {crpix_y} --radius 5" \
              f"  --cpulimit 600 --no-plots --extension 1 {frame_path} "
# --sigma 0.5,  --no-tweak
    subprocess.call(cmd, shell=True)

    wcs_path = f"/home/bekah/glcat/astrometry/aspect_correction/frame{i}.wcs"

    return wcs_path


def get_stars(image, frame, threshold: float = 0.75, star_size: int = 2):
    """ run DAO on movie frames, sort by flux """
    from photutils import DAOStarFinder

    daofind = DAOStarFinder(fwhm=star_size, threshold=threshold, sharplo=0.00)
    star_list = daofind(image)
    if star_list is not None:
        print(f"{len(star_list)} stars found")
        table_name = make_fits_table(star_list.to_pandas(), frame)
        return table_name
    else:
        print("no stars found for this frame.")
        return None


def run_astrometry_net(eclipse, expt, xylist_path, image_width, image_height, ra, dec, crpix_x, crpix_y):
    """ Main call to astrometry.net:
    solve-field: command to solve for aspect solution
    --overwrite: write over files
    --no-tweak: no SIPs correction
    --no-plots: hypothetically, no plots
    -w image width, -e image height
    -L and -H: pixel scale upper and lower bounds, -u app is unit
    -3 and -4: ra and dec
    --crpix-x and --crpix-y: center pixel of image for ra and dec
    --radius: area to search around ra and dec in degrees
    --verify: if there is an existing wcs, verify it
    --verify-ext 1: hdu extension for wcs existing
    then at the end is the path to the x,y list of sources in a FITS table """

    print("running astrometry.net on individual frame of movie")
    print(f"frame is near {ra}, {dec}.")
    cmd = f"solve-field --overwrite --no-plots -w {image_width} -e {image_height} -L 1.2 -H 1.8 -u app" \
          f" -3 {ra} -4 {dec} --crpix-x {crpix_x} --crpix-y {crpix_y} --radius 3" \
          f" --verify '/home/bekah/glcat/astrometry/e00{eclipse}/e{eclipse}-nd-{expt}s-0-f0000-rice.fits' --verify-ext 1 " \
          f"'{xylist_path}'"
    # --no-tweak

    subprocess.call(cmd, shell=True)

    return None


def zero_flag_and_edge(cnt, flag, edge):
    cnt[~np.isfinite(cnt)] = 0
    cnt[np.nonzero(flag)] = 0
    cnt[np.nonzero(edge)] = 0
    return cnt


def make_fits_table(star_list, frame):
    print("making fits table of xy positions")
    star_list = star_list.sort_values(by="flux", ascending=False)
    colx = fits.Column(name='X', format='E', array=star_list['xcentroid'])
    coly = fits.Column(name='Y', format='E', array=star_list['ycentroid'])
    hdu = fits.BinTableHDU.from_columns([colx, coly])
    tableName = f'/home/bekah/glcat/astrometry/aspect_correction/frame{frame}.xyls'
    hdu.writeto(tableName, overwrite=True)
    return tableName


def get_ra_dec(eclipse):
    """ Loading aspect table, set the ra and dec to 0 in the pipeline to get the movie,
     so we have to work around that here to get an initial ra / dec guess for
      astrometry.net. """
    parq = parquet.read_table('/home/bekah/gphoton_working/gPhoton/aspect/aspect.parquet')
    aspect = parq.to_pandas()
    ra = np.mean(aspect[aspect["eclipse"] == eclipse]["ra"])
    dec = np.mean(aspect[aspect["eclipse"] == eclipse]["dec"])

    return ra, dec


def get_aspect_from_wcs(wcs_path):

    frame_wcs = fits.open(wcs_path)
    new_center_ra = frame_wcs[0].header["CRVAL1"]
    new_center_dec = frame_wcs[0].header["CRVAL2"]

    return new_center_ra, new_center_dec
