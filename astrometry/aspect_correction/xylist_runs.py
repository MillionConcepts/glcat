# functions to run astrometry.net with xylists of stars


import pandas as pd
from astropy.io import fits
import subprocess
import os.path
from util import zero_flag_and_edge, get_aspect_from_wcs, get_ra_dec
from image_runs import run_image_frame


def run_xylist(eclipse, expt, num_frames, file_names, dose, **opt):
    """Handler for running astrometry.net with xylists, our main way to run it.
    This means it only intakes picture info, as well as a list of star x, y
    positions."""
    aspect_solution = {}
    failed_frames = []

    for i in range(num_frames):
        print(f"opening frame {i} fits file")
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
          f"'{xylist_path}'"
    # --no-tweak
    # f" --verify '/home/bekah/glcat/astrometry/e{eclipse}/e{eclipse}-nd-{expt}s-0-f0000-rice.fits' --verify-ext 1 " \
    # -L 1.2 -H 1.8 -u app

    subprocess.call(cmd, shell=True)

    return None


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


def make_fits_table(star_list, frame, file_names):
    """ writes fits table of source locations, sorted by flux (highest flux = first)"""
    #TODO: change to return something else
    print("making fits table of xy positions")
    star_list = star_list.sort_values(by="flux", ascending=False)
    colx = fits.Column(name='X', format='E', array=star_list['xcentroid'])
    coly = fits.Column(name='Y', format='E', array=star_list['ycentroid'])
    hdu = fits.BinTableHDU.from_columns([colx, coly])
    tableName = file_names['xylist'][frame]
    hdu.writeto(tableName, overwrite=True)
    return tableName
