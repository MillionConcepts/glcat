"""functions to run astrometry.net with xylists of stars"""

from astropy.io import fits
import subprocess
from aspect_correction.util import get_aspect_from_wcsinfo
import os.path
from hostess.subutils import Viewer


def refine_normal_frame(frame_series, xylist):
    """ main refined for a "normal" frame. extracts stars via daofind, sends
    xylist to astrometry.net, uses wcsinfo to get resulting aspect soln """
    if xylist:
        print(f"extracting sources from frame")
        get_stars(frame_series)
        if os.path.isfile(frame_series["xylist_path"]):
            image_width, image_height, crpix_x, crpix_y = get_image_size(frame_series)
            done = astrometry_hostess_run(
                frame_series["xylist_path"],
                frame_series["aspect_output"],
                image_width,
                image_height,
                frame_series['ra'],
                frame_series['dec'],
                crpix_x,
                crpix_y)
        else:
            print("xylist file doesn't exist.")
            return None
    if not xylist:
        if os.path.isfile(frame_series["backplane_path"]):
            backplane = frame_series['backplane_path']
            done = astrometry_hostess_image_run(
                backplane,
                frame_series["aspect_output"],
                frame_series['ra'],
                frame_series['dec'])
        else:
            print(f"backplane file, "
                  f"{frame_series['backplane_path']} doesn't exist.")
            return None
    print("Getting aspect from wcs file")
    if done:
        try:
            aspect = get_aspect_from_wcsinfo(frame_series['wcs_path'])
            print("got aspect from wcsinfo.")
            return aspect
        finally:
            return None
    return

def get_stars(frame_series):
    """ run DAO on movie frames, sort by flux """
    from photutils import DAOStarFinder
    print(f"Opening frame image: {frame_series['backplane_path']}.")
    image = fits.open(frame_series['backplane_path'])
    daofind = DAOStarFinder(fwhm=2, threshold=0.65, sharplo=0.00)
    star_list = daofind(image[0].data)
    if star_list is not None:
        print(f"{len(star_list)} stars found")
        make_fits_table(star_list.to_pandas(), frame_series['xylist_path'])
        return
    else:
        print("no stars found for this frame.")
        return None


def make_fits_table(star_list, table_name):
    """ writes fits table of source locations, sorted by flux (highest flux = first)"""
    print("making fits table of xy positions")
    star_list = star_list.sort_values(by="flux", ascending=False)
    colx = fits.Column(name='X', format='E', array=star_list['xcentroid'])
    coly = fits.Column(name='Y', format='E', array=star_list['ycentroid'])
    hdu = fits.BinTableHDU.from_columns([colx, coly])
    hdu.writeto(table_name, overwrite=True)
    return


def get_image_size(frame_series):
    """ use fits header to get width / height """
    img = fits.open(frame_series['backplane_path'])
    image_width = img[0].header['NAXIS1']
    image_height = img[0].header['NAXIS2']
    crpix_x = image_width/2  # RA -- TAN
    crpix_y = image_height/2  # DEC -- TAN
    return image_width, image_height, crpix_x, crpix_y


def astrometry_hostess_run(
        xylist_path,
        output_path,
        image_width,
        image_height,
        ra,
        dec,
        crpix_x,
        crpix_y):
    """ use hostess to managae runs to astrometry. keeps things
    formatted nicely. astrometry is done running when solve_process.wait()
    finishes and solve_process.done returns True (hopefully). """
    solve_process = Viewer.from_command(
        "solve-field",
        xylist_path,
        overwrite=True,
        no_plots=True,
        dir_=output_path,
        width=image_width,
        height=image_height,
        scale_units="arcsecperpix",
        scale_low=1.0,
        scale_high=1.6,
        N="none",
        U="none",
        B="none",
        M="none",
        R="none",
        ra=ra,
        dec=dec,
        radius=5,
        temp_axy=True,
        crpix_x=crpix_x,
        crpix_y=crpix_y)
    solve_process.wait()
    return solve_process.done


def astrometry_hostess_image_run(
        backplanepath,
        output_path,
        ra,
        dec):
    """ use hostess to managae runs to astrometry with images
     not xylists """
    print(f"running astrometry on image {backplanepath}")
    solve_process = Viewer.from_command(
        "solve-field",
        backplanepath,
        overwrite=True,
        no_plots=True,
        dir_=output_path,
        scale_units="arcsecperpix",
        scale_low=1.0,
        scale_high=1.6,
        N="none",
        U="none",
        B="none",
        M="none",
        R="none",
        ra=ra,
        dec=dec,
        radius=5,
        temp_axy=True)
    solve_process.wait()
    return solve_process.done
