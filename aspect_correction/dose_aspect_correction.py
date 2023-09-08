"""functions to run astrometry.net with xylists of stars"""

from astropy.io import fits
import subprocess
from util import zero_flag_and_edge, get_aspect_from_wcsinfo
import os.path


def refine_normal_frame(frame_series):
    print(f"extracting sources from frame")
    get_stars(frame_series)
    if os.path.isfile(frame_series["xylist_path"]):
        image_width, image_height, crpix_x, crpix_y = get_image_size(frame_series)
        run_astrometry_net(frame_series['eclipse'], frame_series["xylist_path"], image_width,
                           image_height, frame_series['ra'], frame_series['dec'], crpix_x, crpix_y)
        aspect = get_aspect_from_wcsinfo(wcs_path)
    else:
        aspect = None
    return aspect


def run_astrometry_net(eclipse, xylist_path, output_path, image_width, image_height, ra,
                       dec, crpix_x, crpix_y):
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

    print(f"Running astrometry.net on frame of {eclipse}.")
    cmd = f"solve-field --overwrite --no-plots --dir {output_path} -w {image_width} -e {image_height} " \
          f"--scale-units arcsecperpix --scale-low 1.48 --scale-high 1.52 " \
          f"-N none -U none --axy none -B none -M none -R none " \
          f" -3 {ra} -4 {dec} --crpix-x {crpix_x} --crpix-y {crpix_y} --radius 5 {xylist_path}"
    # --no-tweak
    # f" --verify '/home/bekah/glcat/astrometry/e{eclipse}/e{eclipse}-nd-{expt}s-0-f0000-rice.fits' --verify-ext 1 " \
    # -L 1.2 -H 1.8 -u app
    # RADIUS USED TO BE 3

    subprocess.call(cmd, shell=True)

    return None


def get_stars(frame_series):
    """ run DAO on movie frames, sort by flux """
    from photutils import DAOStarFinder
    image = fits.open(frame_series['backplane_path'])
    daofind = DAOStarFinder(fwhm=2, threshold=0.75, sharplo=0.00)
    star_list = daofind(image)
    if star_list is not None:
        print(f"{len(star_list)} stars found")
        make_fits_table(star_list.to_pandas(), frame_series['xylist_path'])
        return
    else:
        print("no stars found for this frame.")
        return None


def make_fits_table(star_list, tableName):
    """ writes fits table of source locations, sorted by flux (highest flux = first)"""
    print("making fits table of xy positions")
    star_list = star_list.sort_values(by="flux", ascending=False)
    colx = fits.Column(name='X', format='E', array=star_list['xcentroid'])
    coly = fits.Column(name='Y', format='E', array=star_list['ycentroid'])
    hdu = fits.BinTableHDU.from_columns([colx, coly])
    hdu.writeto(tableName, overwrite=True)
    return


def get_image_size(frame_series):
    """ use fits header to get width / height """
    img = fits.open(frame_series['backplane_path'])
    image_width = img[0].header['NAXIS1']
    image_height = img[0].header['NAXIS2']
    crpix_x = image_width/2  # RA -- TAN
    crpix_y = image_height/2  # DEC -- TAN
    return image_width, image_height, crpix_x, crpix_y