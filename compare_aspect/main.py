
import subprocess
import os
import sys
from typing import Optional
import aspect_correction as asp
from backplanes import make_backplanes
from compare_aspect.psf_fitting import run_psf_compare
from compare_aspect.plots import make_plots
sys.path.append('/home/bekah/gphoton_working')
sys.path.append('/home/bekah/glcat/aspect_correction')
sys.path.append('/home/bekah/gphoton_working/gPhoton')


def run_compare(
        eclipse,
        band,
        runtype: Optional[str] = None,
        runnote=""
):
    """run comparison of a given aspect soln (usually the og pipeline one)
     and new aspect solution type, takes as an argument the type
    of aspect improvement to use (astrometry.net, astroalign, etc)
    Args:
     eclipse:eclipse number
     band: "NUV" or "FUV"
     runtype: None or with "short" you can specify if you only
     want to run aspect refiner and psf fitting on the resulting
     image using a pre-existing xy / photom list (just sidesteps
     first gphoton run & psf fitting), or "psf_only" only does
     psf fitting on existing files (must be properly named)
     runnote: optional string to add to the end of output filenames
     """
    print(f"Starting comparison run for eclipse {eclipse}, {band}.")
    file_names = make_file_names(eclipse, runnote)
    # run gphoton with old aspect solution parquet file
    # save output files to "test_datadda" and sub eclipse folder
    if runtype != "psf_only":
        if runtype != "short":
            print(f"Running gphoton with old aspect solution file.")
            run_gphoton(
                eclipse,
                band,
                "/home/bekah/gphoton_working/test_data/",
                "aspect"
            )
        exptime = get_expt(file_names)
        # run backplanes to produce dosemaps or xylists
        # currently just xylists to save space
        print("writing dosemap backplanes")
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
            star_size=2
        )
        # produce new aspect solution
        print("running aspect refiner")
        asp.main.execute_refiner(eclipse, 1, exptime, 'xylist', dose=True, crop=True)
        # write new aspect solution to aspect2.parquet file
        print("writing aspect to aspect parquet 2")
        write_aspect2(eclipse, file_names)
        # run gphoton with new aspect solution parquet file (aspect2)
        # save output files to "astrom_test_data" and sub eclipse folder
        print(f"Running gphoton with new aspect solution file.")
        run_gphoton(
            eclipse,
            band,
            "/home/bekah/gphoton_working/astrom_test_data/",
            "aspect2"
        )
        # check outputs to make sure it worked
        if check_outputs(file_names):
            return f"gphoton failed to produce the expected images."
    # produce plots comparing images and aspect solutions produced by both runs
    make_plots(file_names)
    # fit PSF and compare FWHM of how ever many stars in full-depth images produced
    # by each run
    print("running psf fitting")
    psf_comparison_tab = run_psf_compare(file_names, runtype)
    print(f"saving psf fitting results to csv at {file_names['psf_comp']}")
    psf_comparison_tab.to_csv(file_names['psf_comp'])
    return


def run_gphoton(eclipse, band, local_root, aspect):
    """use gphoton's cli hooks to run and save to a designated directory w/
    unique names"""
    cmd = f"python /home/bekah/gphoton_working/pipeline_cli.py {eclipse} {band}" \
          f" --threads=4 --local_root={local_root}" \
          f" --compression=rice --extended_photonlist=True --aspect={aspect} --verbose=4"
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
    """takes refined aspect table df and turns it into an aspect2 parquet
     file formatted the same as the og aspect parquet file so it can be
      read by gphoton, overwrites whatever is previously called
      aspect2.parquet in the gphoton aspect folder"""
    import pandas as pd
    # import pyarrow
    # check if refined aspect table exists
    if os.path.exists(file_names['new_aspect']):
        new_aspect_df = pd.read_csv(file_names['new_aspect'])
        df2 = new_aspect_df.filter(['ra_center',
                                    'dec_center',
                                    'orientation'],
                                   axis=1)
        df2['eclipse'] = eclipse
        # have to rename columns to ra, dec, and roll
        # bc that's how og aspect parq is
        df2 = df2.rename(columns={"ra_center": "ra",
                                  "dec_center": "dec",
                                  "orientation": "roll"})
        # save as parquet
        df2.to_parquet(file_names["aspect_parq"], compression=None)
    elif not os.path.exists(file_names['new_aspect']):
        print("new aspect solution does not exist.")
    return


def make_file_names(eclipse, runnote=""):
    """make file names for running aspect comparison. will need to be updated
    should not-Bekah need to use, or if gphoton file names change bc of diff
    run params.
    Args:
        eclipse: the eclipse number of the run, will get padded
        runnote: optional way to add text to output file names
         to distinguish diff kinds of runs on the fly """
    padded_eclipse = str(eclipse).zfill(5)
    # set up dirs and make dirs if they don't exist
    base = "/home/bekah/gphoton_working/"
    oldbase = f"{base}test_data/e{padded_eclipse}/"
    newbase = f"{base}astrom_test_data/e{padded_eclipse}/"
    compbase = f"{base}test_data_comp/e{padded_eclipse}/"
    if not os.path.exists(compbase):
         os.mkdir(compbase)
    # outputs images
    old_image_file = oldbase+f"e{padded_eclipse}" \
                             f"-nd-tfull-b00-image-r.fits"
    new_image_file = newbase+f"e{padded_eclipse}" \
                             f"-nd-tfull-b00-image-r.fits"
    image_comparison = compbase+f"e{padded_eclipse}" \
                       f"+image_comp"+runnote+".jpg"
    old_photom_file = oldbase+f"e{padded_eclipse}" \
                      f"-nd-tfull-b00-image-photom-12_8.csv"
    new_aspect = newbase+\
                 f"astrometry_xy_1s_{padded_eclipse}"
    ra_dec_plot = compbase+\
                  f"e{padded_eclipse}_ra_dec_comp"+runnote+".jpg"
    old_asp_title = f"e{eclipse} old asp"
    new_asp_title = f"e{eclipse} new asp"
    star_cutouts = compbase+f"e{eclipse}_star_cutouts"+runnote+".jpg"
    psf_comp = compbase+f"e{eclipse}_psf"+runnote+".csv"
    aspect_parq = base + "gPhoton/aspect/aspect2.parquet"
    file_names = {'old_image_file': old_image_file,
                  'new_image_file': new_image_file,
                  'image_comparison': image_comparison,
                  'old_photom_file': old_photom_file,
                  'new_aspect': new_aspect,
                  'ra_dec_plot': ra_dec_plot,
                  'star_cutouts': star_cutouts,
                  'old_asp_title': old_asp_title,
                  'new_asp_title': new_asp_title,
                  'psf_comp': psf_comp,
                  'aspect_parq': aspect_parq}
    return file_names


def get_expt(file_names):
    """opens fits image file and looks at header to get expt"""
    from astropy.io import fits
    image = fits.open(file_names['old_image_file'])
    expt = image[1].header["EXPTIME"]
    return int(expt)


