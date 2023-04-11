
from compare_aspect.plots import make_plots
import subprocess
import os
import astrometry.aspect_correction as asp
from backplanes import make_backplanes
import sys

sys.path.append('/home/bekah/gphoton_working')
sys.path.append('/home/bekah/glcat/astrometry/aspect_correction')
sys.path.append('/home/bekah/gphoton_working/gPhoton')


def run_compare(eclipse, band, old_aspect, new_aspect_method):
    """run comparison of a given aspect soln (usually the og pipeline one)
     and new aspect solution type, takes as an argument the type
    of aspect improvement to use (astrometry.net, astroalign, etc)
    Args:
     eclipse:eclipse number
     band: "NUV" or "FUV"
     old_aspect:o.g. pipeline, could point to a diff method?
     new_aspect_method:'astrometry', 'astroalign', 'RANSAC', 'drizzlepac'
     """
    print(f"Starting comparison run for eclipse {eclipse}, {band}.")
    file_names = make_file_names(eclipse, old_aspect, new_aspect_method)
    # run gphoton with old aspect solution parquet file
    # save output files to "test_data" and sub eclipse folder
    print(f"Running gphoton with old aspect solution file.")
    run_gphoton(eclipse, band, "test_data", "aspect")
    # run backplanes to produce dosemaps or xylists
    # currently just xylists to save space
    make_backplanes(
        eclipse=eclipse,
        band=band,
        depth=1,
        leg=0,
        threads=4,
        burst=True,
        local="/home/bekah/gphoton_working/test_data",
        kind="dose",
        radius=400,
        write={'array': False, 'xylist': True},
        inline=True,
        threshold=.75,
        star_size=2)
    # produce new aspect solution
    asp.main.execute_refiner(eclipse, 1, 150, 'xylist', dose=True, crop=False)
    # write new aspect solution to aspect2.parquet file
    write_aspect2(eclipse, file_names)
    # run gphoton with new aspect solution parquet file (aspect2)
    # save output files to "astrom_test_data" and sub eclipse folder
    print(f"Running gphoton with new aspect solution file.")
    run_gphoton(eclipse, band, "astrom_test_data", "aspect2")
    # check outputs to make sure it worked
    if check_outputs(file_names):
        return f"gphoton failed to produce the expected images."
    # produce plots comparing images and aspect solutions produced by both runs
    make_plots(eclipse, old_aspect, new_aspect_method, file_names)
    # fit PSF and compare FWHM of how ever many stars in full-depth images produced
    # by each run
    return


def run_gphoton(eclipse, band, local_root, aspect):
    """use gphoton's cli hooks to run and save to a designated directory w/
    unique names"""
    cmd = f"python pipeline_cli.py {eclipse} {band} --threads=4 --local_root={local_root}" \
          f"--compression=rice --extended_photonlist=True --aspect={aspect}"
    subprocess.call(cmd, shell=True)

    return


def check_outputs(file_names):
    """check that gphoton produced the expected files"""
    if os.path.exists(file_names['old_image_file']):
        if os.path.exists(file_names['new_image_file']):
            print("Both full-depth images created successfully in gPhoton.")
            return False
    print("One or both full-depth images failed in gPhoton.")
    return True


def write_aspect2(eclipse, file_names):
    """takes refined aspect table df and turns it into an aspect2 parquet file
     formatted the same as the og aspect parquet file so it can be read by gphoton,
      overwrites whatever is previously called aspect2.parquet in the gphoton aspect
      folder"""
    import pandas as pd
    # import pyarrow
    # check if refined aspect table exists
    if os.path.exists(file_names['new_aspect']):
        new_aspect_df = pd.read_csv(file_names['new_aspect'])
        df2 = new_aspect_df.filter(['ra_center', 'dec_center', 'orientation'], axis=1)
        df2['eclipse'] = eclipse
        # have to rename columns to ra, dec, and roll bc that's how og aspect parq is
        df2.rename(columns={"ra_center": "ra", "dec_center": "dec", "orientation": "roll"})
        # save as parquet
        df2.to_parquet('/home/bekah/gphoton_working/gPhoton/aspect2.parquet', compression=None)


def make_file_names(eclipse, old_aspect, new_aspect_method):
    oldbase = f"/home/bekah/gphoton_working/test_data/e{eclipse}/"
    newbase = f"/home/bekah/gphoton_working/astrom_test_data/e{eclipse}/"
    compbase = f"/home/bekah/gphoton_working/test_data_comp/e{eclipse}/"
    padded_eclipse = str(eclipse).zfill(5)
    old_image_file = oldbase+f"e{padded_eclipse}-nd-tfull-b00-image-r.fits"
    new_image_file = newbase+f"e{padded_eclipse}-nd-tfull-b00-image-r.fits"
    image_comparison = compbase+f"e{eclipse}+image_comp.jpg"
    old_aspect = oldbase
    new_aspect = newbase+f"astrometry_xy_1s_{eclipse}"
    ra_dec_plot = compbase+f"e{eclipse}_ra_dec_comp.jpg"
    old_asp_title = f"e{eclipse} old asp"
    new_asp_title = f"e{eclipse} new asp"
    file_names = {'old_image_file': old_image_file, 'new_image_file': new_image_file,
                  'image_comparison': image_comparison, 'old_aspect': old_aspect,
                  'new_aspect': new_aspect, 'ra_dec_plot': ra_dec_plot,
                  'old_asp_title': old_asp_title, 'new_asp_title': new_asp_title}
    return file_names



