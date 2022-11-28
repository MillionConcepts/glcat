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


def execute_refiner(eclipse):
    print("opening fits file")
    movie = fits.open(f"/home/bekah/glcat/astrometry/e23456/e{eclipse}-nd-10s-rice.fits")

   # cnt, flag, edge = [hdu.read() for hdu in movie[1:4]]

    print("masking edge and flags")
    #masked_cnt_image = zero_flag_and_edge(cnt, flag, edge)

    num_frames = len(movie[1].data)
    print(f"there are {num_frames} frames in the movie.")

    aspect_solution = {}
    for i in range(num_frames):
        print(f"extracting sources from frame {i}")
        table_name = get_stars(movie[1].data[i], i)
        if table_name is None:
            print(f"No stars found in this frame {i}.")
        else:
          image_width = movie[1].header['NAXIS1']
          image_height = movie[1].header['NAXIS2']
          crpix_x = movie[1].header['CRPIX1'] # RA -- TAN
          crpix_y = movie[1].header['CRPIX2'] # DEC -- TAN

          ra, dec = get_ra_dec(eclipse)

          run_astrometry_net(eclipse, table_name, image_width, image_height, ra, dec, crpix_x, crpix_y)

          new_center_ra, new_center_dec = get_aspect_from_wcs(i)

          aspect_solution[i] = [new_center_ra, new_center_dec]

    asp_df = pd.DataFrame(aspect_solution)
    asp_df.to_csv(f"{eclipse}_aspect_soln")

    return aspect_solution

def get_stars(image, frame):
    """ run DAO on movie frames, sort by flux """
    from photutils import DAOStarFinder

    threshold = 1
    daofind = DAOStarFinder(fwhm=3, threshold=threshold, sharplo=0.00)
    star_list = daofind(image)
    if star_list is not None:
        table_name = make_fits_table(star_list.to_pandas(), frame)
        return table_name
    else:
        return None

def run_astrometry_net(eclipse, xylist_path, image_width, image_height, ra, dec, crpix_x, crpix_y):
    print("running astrometry.net on individual frame of movie")
    print(f"frame is near {ra}, {dec}.")
    cmd = f"solve-field --overwrite --no-plots -w {image_width} -e {image_height} -L 1.2 -H 1.6 -u app" \
          f" -3 {ra} -4 {dec} --crpix-x crpix_x --crpix-y crpix_y --radius 3" \
          f" --verify '/home/bekah/glcat/astrometry/e23456/e{eclipse}-nd-10s-rice.fits' --verify-ext 1 " \
          f"'{xylist_path}'"

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
    coldefs = fits.ColDefs([colx, coly])
    hdu = fits.BinTableHDU.from_columns([colx, coly])
    tableName = f'/home/bekah/aspect_correction/frame{frame}.xyls'
    hdu.writeto(tableName, overwrite=True)
    return tableName

def get_ra_dec(eclipse):
    # loading aspect table
    # we set the ra and dec to 0 in the pipeline to get the movie, so we have to work around that here
    # (I think)
    parq = parquet.read_table('/home/bekah/gphoton_working/gPhoton/aspect/aspect.parquet')
    aspect = parq.to_pandas()
    ra = np.mean(aspect[aspect["eclipse"] == eclipse]["ra"])
    dec = np.mean(aspect[aspect["eclipse"] == eclipse]["dec"])

    return ra, dec

def get_aspect_from_wcs(i):
    # could not do this is a subprocess call by rewriting some of the astrometry code from c, or learning c lol
    cmd = f"frame{i}.wcs"
    aspect = subprocess.run(['wcsinfo', cmd], stdout=subprocess.PIPE).stdout.decode('utf-8').split()
    new_center_ra = aspect[37] # ideally turn aspect into a dict so this is less hacked together
    new_center_dec = aspect[39]
    return new_center_ra, new_center_dec