""" meant to be run after producing a new aspect solution that's a csv. creates a kind of
 'quality' report on the aspect soln """

import os
import sys
from psf_fitting_new import run_psf_compare
from plots import make_plots
import pandas as pd
sys.path.append('/aspect_correction')
sys.path.append('/gPhoton')
sys.path.append('/home/bekah/gPhoton2')
sys.path.append('/home/ubuntu/glcat/aspect_correction')
sys.path.append('/home/ubuntu/gPhoton2/gPhoton')
from gPhoton.pipeline import execute_pipeline


# aspect solution has been produced, is a csv
def run_compare(
        eclipse,
        leg,
        band,
        local_root="/home/bekah/gPhoton2",
        runnote=""
):
    """run comparison of a given aspect soln (usually the og pipeline one)
     and new aspect solution type, takes as an argument the type
    of aspect improvement to use (astrometry.net, astroalign, etc)
    Args:
     eclipse:eclipse number
     band: "NUV" or "FUV"
     local_root: where files go
     runnote: optional string to add to the end of output filenames
     """
    print("making local folders for saving outputs.")
    og_folder, new_folder = make_local_directories(local_root, eclipse)
    file_names = get_file_names(eclipse, band, leg, og_folder, new_folder, local_root)
    print("writing aspect to aspect parquet 2")
    write_aspect2(eclipse, file_names)
    print(f"Starting comparison run for eclipse {eclipse}, {band}.")
    execute_pipeline(
        eclipse,
        band,
        depth=None,
        threads=4,
        local_root=og_folder,
        recreate=True,
        aperture_sizes=[12.8],
        write={"movie": False, "image": True},
        compression="rice",
        lil=True,
        burst=False,
        extended_photonlist=True,
        aspect="aspect"
    )
    print(f"Running gphoton with new aspect solution file for {eclipse}.")
    execute_pipeline(
        eclipse,
        band,
        depth=None,
        threads=4,
        local_root=new_folder,
        recreate=True,
        aperture_sizes=[12.8],
        write={"movie": False, "image": True},
        compression="rice",
        lil=True,
        burst=False,
        extended_photonlist=True,
        aspect="aspect2"
    )
    if check_outputs(file_names):
        return f"gphoton failed to produce the expected images."
    make_plots(file_names)
    print("running psf fitting")
    psf_comparison_tab = run_psf_compare(file_names)
    print(f"saving psf fitting results to csv at {file_names['psf_comp']}")
    psf_comparison_tab.to_csv(file_names['psf_comp'])
    return


def get_file_names(eclipse, band, leg, og_folder, new_folder, root):
    """ image file names and photonlists """
    file_names = {}
    eclipse_num = str(eclipse).zfill(5)
    eclipse_str = "e" + eclipse_num
    leg_str = str(leg).zfill(2)
    if band == "NUV":
        b = "nd"
    elif band == "FUV":
        b = "fd"
    og_folder = og_folder + "/" + eclipse_str + "/"
    new_folder = new_folder + "/" + eclipse_str + "/"
    file_names['old_image_file'] = og_folder + eclipse_str + f"-{b}-b{leg_str}-ffull-image-r.fits"
    file_names['new_image_file'] = new_folder + eclipse_str + f"-{b}-b{leg_str}-ffull-image-r.fits"
    # if aperture changes, 12_8 will change. could propagate up as a variable.
    file_names['psf_comp'] = og_folder + eclipse_str + f"psf_comp.csv"
    file_names['old_photom_file'] = og_folder + eclipse_str + f"-{b}-b{leg_str}-" \
                                                              f"ffull-image-" \
                                                              f"photom-12_8.csv"
    file_names['image_comparison'] = og_folder + eclipse_str+"-image-compare.jpg"
    file_names['new_aspect'] = root+"test_data/"+eclipse_num + "_new_aspect.csv"
    file_names['old_aspect_parq'] = root+'gPhoton/aspect/aspect.parquet'
    file_names['aspect_parq'] = root+'gPhoton/aspect/aspect2.parquet'
    file_names['star_cutouts'] = og_folder + eclipse_str +"-star-cutouts.jpg"
    file_names['star_cutouts_new'] = og_folder + eclipse_str +"-star-cutouts-newasp.jpg"

    return file_names


def check_outputs(file_names):
    """check that gphoton produced the expected files"""
    print(f"Searching for images {file_names['old_image_file']} "
          f"and {file_names['new_image_file']}.")
    if os.path.exists(file_names['old_image_file']):
        if os.path.exists(file_names['new_image_file']):
            return False
    print("One or both full-depth images failed in gPhoton.")
    return True


def write_aspect2(eclipse, file_names):
    """takes refined aspect table df and turns it into an aspect2 parquet
     file formatted the same as the og aspect parquet file so it can be
      read by gphoton, overwrites whatever is previously called
      aspect2.parquet in the gphoton aspect folder. adds time and flag
      info from og aspect."""
    from pyarrow import parquet
    # loading aspect table to add time stamp and flags
    if os.path.exists(file_names['new_aspect']):
        new_aspect_df = pd.read_csv(file_names['new_aspect'])
        df2 = new_aspect_df.filter(['mission_time',
                                    'ra_center',
                                    'dec_center',
                                    'orientation',
                                    'flags_x',
                                    'logodds',
                                    'frame_type',
                                    'time'],
                                   axis=1)
        df2['eclipse'] = eclipse
        # have to rename columns to ra, dec, and roll
        # bc that's how og aspect parq is
        df2 = df2.rename(columns={"ra_center": "ra",
                                  "dec_center": "dec",
                                  "orientation": "roll",
                                  "time": "frame",
                                  "mission_time": "time",
                                  "flags_x": "flags"})
        df2 = df2.set_index('frame')
        df2 = df2.reindex(range(len(df2)), fill_value=0)
        df2 = df2.astype({'ra': 'float64',
                          'dec': 'float64',
                          'roll': 'float64',
                          'time': 'float64',
                          'logodds': 'str',
                          'flags': 'float64',
                          'frame_type': 'str'})
        print("Interpolating df")
        #df2 = df2.interpolate(method='linear',
        #                      limit_direction='both',
        #                      axis=0)
        # save to parquet
        df2.to_parquet(file_names["aspect_parq"], compression=None)
    elif not os.path.exists(file_names['new_aspect']):
        print("new aspect solution does not exist.")
        return
    return


def make_local_directories(local_root, eclipse):
    import os
    og_folder = local_root + "og_asp_"+str(eclipse)
    new_folder = local_root + "new_asp"+str(eclipse)
    if not os.path.exists(og_folder):
        os.makedirs(og_folder)
    if not os.path.exists(new_folder):
        os.makedirs(new_folder)
    print(f"{og_folder} and {new_folder} are now created.")

    return og_folder, new_folder


def image_plot(file_names):
    """ old and new aspect images  """
    from astropy.io import fits
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(2)
    old_image = fits.open(file_names['old_image_file'])
    new_image = fits.open(file_names['new_image_file'])
    axs[0].imshow(centile_clip(old_image[1].data, (0, 99.5)),
                  interpolation='none', cmap='Greys_r', origin='lower')
    axs[1].imshow(centile_clip(new_image[1].data, (0, 99.5)),
                  interpolation='none', cmap='Greys_r', origin='lower')
    plt.savefig(file_names['image_comparison'])
    return


def centile_clip(image, centiles=(0, 90)):
    """
    simple clipping function that clips values above and below a given
    percentile range
    """
    import numpy as np
    finite = np.ma.masked_invalid(image)
    bounds = np.percentile(finite[~finite.mask].data, centiles)
    result = np.ma.clip(finite, *bounds)
    if isinstance(image, np.ma.MaskedArray):
        return result
    return result.data