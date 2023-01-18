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
# 3) Change the output location of automatic astrometry.net files like wcs etc (also why does it make so many :( )
# 4) delete all astrometry.net files after each frame
# 5) try to plot xylists for each frame to see how they evolve (ie are they always in the same spot on the star?)
#    they're probably not
# 6) get rid of random 0 and 1 at the beginning of the refined aspect solution output csv

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
    Main pipeline for aspect refinement on individual frames (short time period, not projected)
    of eclipses. Sources are detected using DAO star finder, source list + intiial ra / dec are fed to
    astrometry.net when using an "xylist" run.
    Inputs:
    * eclipse = for calling files
    * expt = number of seconds per frame used in gphoton2, for file name
    * num_frames = can come up with a better way to do this,
     currently just for knowing how many times to run loop2
    * run_type = can be "image" or "xylist". xylist runs DAOstarfinder
     and feeds list of sources to astrometry.net, image feeds the image
     directly to astrometry.net. image tends to take longer.
    * threshold = threshold for star detection used by dao, optional
    * star_size = what size of star dao is looking for, will probably
      vary with the frame length used in the run
    * verify = runs a second run of astrometry.net for "verification"
     that uses outputs of previous run. may take a while, but hopefully
     improves output.
     """
    print(f"Executing refiner for eclipse {eclipse} using astrometry.net")
    print(f"there are {num_frames} frames in the movie.")

    # setting options and padding eclipse num
    opt = {
        'threshold': threshold,
        'star_size': star_size,
    }
    #eclipse = str(eclipse).zfill(5)

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
        #old_wcs = f'/home/bekah/glcat/astrometry/aspect_correction/frame{i}.wcs'
        if os.path.exists(old_wcs):
            cmd = f"solve-field -v --overwrite --verify {old_wcs} -y -3 {ra} -4 {dec} " \
                  f"-w {image_width} -e {image_height} -L 1.2 -H 1.8 -u app " \
                  f"  --crpix-x {crpix_x} --crpix-y {crpix_y} --radius 5" \
                  f" --cpulimit 500 --extension 1 {frame_path} "
            # --sigma 0.5,  --no-tweak
            subprocess.call(cmd, shell=True)

        #TODO: fix new wcs path ot be in dictionary
        wcs_path = file_names["ver_wcs"][frame] #f"/home/bekah/glcat/astrometry/aspect_correction/e21442/e21442-nd-1s-f{frame_padded}-rice.wcs"

        if os.path.exists(wcs_path):
            new_center_ra, new_center_dec = get_aspect_from_wcs(wcs_path)
            aspect_solution[frame] = [new_center_ra, new_center_dec]

    asp_df = pd.DataFrame(aspect_solution)

    return asp_df


def run_xylist(eclipse, expt, num_frames, file_names, dose, **opt):

    path = f'/home/bekah/glcat/astrometry/e{eclipse}/'

    aspect_solution = {}
    failed_frames = []

    for i in range(num_frames):
        print(f"opening frame {i} fits file")
        #frame_padded = str(i).zfill(4)
        #frame_path = path+f"e{eclipse}-nd-{expt}s-0-f{frame_padded}-rice.fits" # -0, leg
        frame_path = file_names["frame_path"][i]
        movie = fits.open(frame_path)

        if not dose:
            print(f"masking flag and edge in frame {i}")
            masked_cnt = zero_flag_and_edge(movie[1].data[0], movie[2].data[0], movie[3].data[0])
            print(f"extracting sources from frame {i}")
            table_name = get_stars(masked_cnt, i, file_names, **opt)

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
                crpix_x = movie[1].header['CRPIX1']   # RA -- TAN (was -1)
                crpix_y = movie[1].header['CRPIX2']  # DEC -- TAN (was -1)
                print(f"crpix_x = {crpix_x}, crpix_y = {crpix_y}")
                print(f"image width = {image_width}, image height = {image_height}")

                run_astrometry_net(eclipse, expt, table_name, image_width, image_height, ra, dec, crpix_x, crpix_y)

        if dose:
            print(f"extracting sources from frame {i}")
            table_name = get_stars(movie[0].data, i, file_names, **opt)

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
                image_width = 3200
                image_height = 3200
                crpix_x = 1600  # RA -- TAN
                crpix_y = 1600  # DEC -- TAN
                print(f"crpix_x = {crpix_x}, crpix_y = {crpix_y}")
                print(f"image width = {image_width}, image height = {image_height}")

                run_astrometry_net(eclipse, expt, table_name, image_width, image_height, ra, dec, crpix_x, crpix_y)

        wcs_path = f"/home/bekah/glcat/astrometry/aspect_correction/frame{i}.wcs"
        if os.path.exists(wcs_path):
            new_center_ra, new_center_dec = get_aspect_from_wcs(wcs_path)
            aspect_solution[i] = [new_center_ra, new_center_dec]
        if not os.path.exists(wcs_path):
            # try an alternative solution without DAO-produced XY lists
            run_image_frame(eclipse, expt, i, file_names)
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
    cmd = f"solve-field --overwrite --no-plots -w {image_width} -e {image_height} " \
          f"--scale-units arcsecperpix --scale-low 1.4 --scale-high 1.6 " \
          f" -3 {ra} -4 {dec} --crpix-x {crpix_x} --crpix-y {crpix_y} --radius 3" \
          f" --verify '/home/bekah/glcat/astrometry/e{eclipse}/e{eclipse}-nd-{expt}s-0-f0000-rice.fits' --verify-ext 1 " \
          f"'{xylist_path}'"
    # --no-tweak
    # -L 1.2 -H 1.8 -u app

    subprocess.call(cmd, shell=True)

    return None


def run_image(eclipse, expt, num_frames, file_names):

    # TODO: need way to filter out images with only 0s
    #  (typically at beginnning of frames) bc it slows everything down

    path = f'/home/bekah/glcat/astrometry/e{eclipse}/'
    aspect_solution = {}
    ra, dec = get_ra_dec(eclipse)

    for i in range(num_frames):

        frame_path = file_names["frame_path"][i]
        #frame_padded = str(i).zfill(4)
        #frame_path = str(path+f"e{eclipse}-nd-{expt}s-0-f{frame_padded}-rice.fits")
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

        wcs_path = file_names["output_wcs"][i] #f"/home/bekah/glcat/astrometry/aspect_correction/frame{i}.wcs"
        if os.path.exists(wcs_path):
            new_center_ra, new_center_dec = get_aspect_from_wcs(wcs_path)
            aspect_solution[i] = [new_center_ra, new_center_dec]

    asp_df = pd.DataFrame(aspect_solution)
    return asp_df


def run_image_frame(eclipse, expt, i, file_names):
    #  for frame that fails w xy-list

    path = f'/home/bekah/glcat/astrometry/e{eclipse}/'
    aspect_solution = {}
    ra, dec = get_ra_dec(eclipse)

    frame_path = file_names["frame_path"][i]
    #frame_padded = str(i).zfill(4)
    #frame_path = str(path+f"e{eclipse}-nd-{expt}s-0-f{frame_padded}-rice.fits")

    print("running astrometry.net on individual frame of movie")
    print(frame_path)

    movie = fits.open(frame_path)

    image_width = movie[1].header['NAXIS1']
    image_height = movie[1].header['NAXIS2']
    crpix_x = movie[1].header['CRPIX1']   # RA -- TAN
    crpix_y = movie[1].header['CRPIX2']   # DEC -- TAN

    cmd = f"solve-field --overwrite --sigma 0.2 -y -3 {ra} -4 {dec} " \
              f"-w {image_width} -e {image_height} --scale-units arcsecperpix --scale-low 1.4 --scale-high 1.6" \
              f"  --crpix-x {crpix_x} --crpix-y {crpix_y} --radius 5" \
              f"  --cpulimit 600 --no-plots --extension 1 {frame_path} "
# --sigma 0.5,  --no-tweak
    subprocess.call(cmd, shell=True)

    wcs_path = file_names["output_wcs"][i]
    #wcs_path = f"/home/bekah/glcat/astrometry/aspect_correction/frame{i}.wcs"

    return wcs_path


def get_stars(image, frame, file_names, threshold: float = 0.75, star_size: int = 2):
    """ run DAO on movie frames, sort by flux """
    from photutils import DAOStarFinder

    daofind = DAOStarFinder(fwhm=star_size, threshold=threshold, sharplo=0.00)
    star_list = daofind(image)
    if star_list is not None:
        print(f"{len(star_list)} stars found")
        table_name = make_fits_table(star_list.to_pandas(), frame, file_names)
        return table_name
    else:
        print("no stars found for this frame.")
        return None


#TODO: Helper functions below, should be moved to a different file eventually


def make_file_names(eclipse, expt, num_frames, run_type, threshold, star_size):
    """ make dict of file names for output files produced in aspect refinement
    (can also be used to delete files later, because we don't need to keep most
     of these). not all of these names will be used for every run time, but
      strings don't take up a lot of room. """
    dose = True

    fits_table = []
    frame_path = []
    output_wcs = []
    ver_wcs = []

    if not dose:
        image_path = f'/home/bekah/glcat/astrometry/e{eclipse}/'
        for i in range(num_frames):
            # for feeding xylist to astrometry.net
            # think there's a fancier list comprehension way to do this but can't google bc on a plane
            fits_table.append(f'/home/bekah/glcat/astrometry/aspect_correction/frame{i}.xyls')
            # paths to movie frame images, largely used for an image or verification run
            frame_padded = str(i).zfill(4)
            frame_path.append(str(image_path + f"e{eclipse}-nd-{expt}s-f{frame_padded}-rice.fits")) # 0-
            # for each frame, resulting new wcs
            output_wcs.append(f"/home/bekah/glcat/astrometry/aspect_correction/frame{i}.wcs")
            # verified wcs list
            ver_wcs.append(f"/home/bekah/glcat/astrometry/aspect_correction/e21442/e21442-nd-1s-f{frame_padded}-rice.wcs") #TODO: don't hardcode this
        # may not be used, for output of a verification run
        verified_df = f"{eclipse}_aspect_soln_{expt}s_thresh{threshold}_size{star_size}_verified"
        # main refined aspect solution
        if threshold and star_size is None:
            asp_df = f"{eclipse}_aspect_soln_{expt}s"
        else:
            asp_df = f"{eclipse}_aspect_soln_{expt}s_thresh{threshold}_size{star_size}_noMinus1"

    if dose:
        image_path = f'/home/bekah/glcat/astrometry/e{eclipse}/dose_files/'
        for i in range(num_frames):
            # for feeding xylist to astrometry.net
            # think there's a fancier list comprehension way to do this but can't google bc on a plane
            fits_table.append(f'/home/bekah/glcat/astrometry/aspect_correction/frame{i}.xyls')
            # paths to movie frame images, largely used for an image or verification run
            frame_padded = str(i).zfill(5)
            time_padded = str(expt).zfill(4)
            frame_path.append(str(image_path + f"e{eclipse}-nd-t{time_padded}-b00-f{frame_padded}-g_dose.fits")) # 0-
            # for each frame, resulting new wcs
            output_wcs.append(f"/home/bekah/glcat/astrometry/aspect_correction/frame{i}.wcs")
            # verified wcs list
            ver_wcs.append(f"/home/bekah/glcat/astrometry/aspect_correction/e21442/e21442-nd-1s-f{frame_padded}-rice.wcs") #TODO: don't hardcode this
        # may not be used, for output of a verification run
        verified_df = f"{eclipse}_aspect_soln_{expt}s_thresh{threshold}_size{star_size}_verified"
        # main refined aspect solution
        if threshold and star_size is None:
            asp_df = f"{eclipse}_aspect_soln_{expt}s"
        else:
            asp_df = f"{eclipse}_aspect_soln_{expt}s_thresh{threshold}_size{star_size}_dose"

    file_names = {"asp_df": asp_df, "verified_df": verified_df, "xylist": fits_table,
                  "frame_path": frame_path, "output_wcs": output_wcs, "ver_wcs": ver_wcs}
    return file_names


def zero_flag_and_edge(cnt, flag, edge):
    """ set flag and edge pixels to zero within count array, used to prevent
    fake sources. """
    cnt[~np.isfinite(cnt)] = 0
    cnt[np.nonzero(flag)] = 0
    cnt[np.nonzero(edge)] = 0
    return cnt


def make_fits_table(star_list, frame, file_names):
    """ writes fits table of source locations, sorted by flux (highest flux = first)"""
    #TODO: change to return something else
    print("making fits table of xy positions")
    star_list = star_list.sort_values(by="flux", ascending=False)
    colx = fits.Column(name='X', format='E', array=star_list['xcentroid'])
    coly = fits.Column(name='Y', format='E', array=star_list['ycentroid'])
    hdu = fits.BinTableHDU.from_columns([colx, coly])
    tableName = file_names['xylist'][frame] #f'/home/bekah/glcat/astrometry/aspect_correction/frame{frame}.xyls'
    hdu.writeto(tableName, overwrite=True)
    return tableName


def get_ra_dec(eclipse):
    """ Loading aspect table, set the ra and dec to 0 in the pipeline to get the movie,
     so we have to work around that here to get an initial ra / dec guess for
      astrometry.net. """
    #eclipse = 9869
    parq = parquet.read_table('/home/bekah/gphoton_working/gPhoton/aspect/aspect.parquet')
    aspect = parq.to_pandas()
    ra = np.mean(aspect[aspect["eclipse"] == eclipse]["ra"])
    dec = np.mean(aspect[aspect["eclipse"] == eclipse]["dec"])
    return ra, dec


def get_aspect_from_wcs(wcs_path):
    """To get the new ra and dec for each frame, open the output wcs file."""
    #TODO: it would be nice to not have to close and open so many files, also check
    # that these files actually each close after use
    frame_wcs = fits.open(wcs_path)
    new_center_ra = frame_wcs[0].header["CRVAL1"]
    new_center_dec = frame_wcs[0].header["CRVAL2"]
    return new_center_ra, new_center_dec

