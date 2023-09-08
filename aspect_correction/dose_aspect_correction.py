"""functions to run astrometry.net with xylists of stars"""

from astropy.io import fits
import subprocess
from util import zero_flag_and_edge, get_aspect_from_wcsinfo


def run_xylist_on_image(eclipse, expt, num_frames, file_names, dose, **opt):
    """Handler for running astrometry.net with xylists, our main way to run it.
    This means it only intakes picture info, as well as a list of star x, y
    positions."""
    for i in range(num_frames):
        print(f"opening frame {i} fits file")
        frame_path = file_names["frame_path"][i]
        movie = fits.open(frame_path)
        ra, dec = get_ra_dec(eclipse)
        if not dose:
            print(f"masking flag and edge in frame {i}")
            masked_cnt = zero_flag_and_edge(movie[1].data[0],
                                            movie[2].data[0],
                                            movie[3].data[0])
            print(f"extracting sources from frame {i}")
            table_name = get_stars(masked_cnt, i, file_names, **opt)
            if table_name is None:
                print(f"No stars found in this frame {i}.")
            else:
                image_width = movie[1].header['NAXIS1']
                image_height = movie[1].header['NAXIS2']
                crpix_x = movie[1].header['CRPIX1']   # RA -- TAN (was -1)
                crpix_y = movie[1].header['CRPIX2']  # DEC -- TAN (was -1)
                run_astrometry_net(eclipse, i, expt, table_name, image_width,
                                   image_height, ra, dec, crpix_x, crpix_y)
        if dose:
            print(f"extracting sources from frame {i}")
            table_name = get_stars(movie[0].data, i, file_names, **opt)
            if table_name is None:
                print(f"No stars found in this frame {i}.")
            if table_name is not None:
                image_width = 3200
                image_height = 3200
                crpix_x = 1600  # RA -- TAN
                crpix_y = 1600  # DEC -- TAN
                run_astrometry_net(eclipse, i, expt, table_name, image_width,
                                   image_height, ra, dec, crpix_x, crpix_y)
    asp_df = make_refined_aspect_table(file_names)
    return asp_df


def run_astrometry_net(eclipse, i, expt, file_names, image_width, image_height, ra,
                       dec, crpix_x, crpix_y, xylist_type):
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

    print(f"Running astrometry.net on {expt} s frame of {eclipse}.")
    print(f"Frame is near {ra}, {dec}.")
    print(f"xylist path is: {file_names[xylist_type][i]}")
    cmd = f"solve-field --overwrite --no-plots --dir {file_names['astrometry_temp']} -w {image_width} -e {image_height} " \
          f"--scale-units arcsecperpix --scale-low 1.48 --scale-high 1.52 " \
          f"-N none -U none --axy none -B none -M none -R none " \
          f" -3 {ra} -4 {dec} --crpix-x {crpix_x} --crpix-y {crpix_y} --radius 5 {file_names[xylist_type][i]}"
    # --no-tweak
    # f" --verify '/home/bekah/glcat/astrometry/e{eclipse}/e{eclipse}-nd-{expt}s-0-f0000-rice.fits' --verify-ext 1 " \
    # -L 1.2 -H 1.8 -u app
    # RADIUS USED TO BE 3

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
