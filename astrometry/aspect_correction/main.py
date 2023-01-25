""" Main pipeline for aspect refinement on individual frames (short time period, not projected)
    of eclipses. Sources are detected using DAO star finder, source list + intiial ra / dec
    are fed to astrometry.net when using an "xylist" run.

    Uses astrometry.net (https://astrometry.net/) takes astronomical images and uses "catalogs"
     to solve for refined aspect solutions. """

from typing import Optional
from xylist_runs import run_xylist_on_image, run_xylist
from image_runs import run_image_frame, run_image, run_verification
from util import make_file_names


def execute_refiner(
                    eclipse,
                    expt,
                    num_frames,
                    run_type,
                    dose: Optional[bool] = False,
                    threshold: Optional[float] = None,
                    star_size: Optional[int] = None,
                    output_directory: Optional[str] = "/home/bekah/glcat_tests/"
                    ):
    """
    Args:
        eclipse: for calling files
        expt: number of seconds per frame used in gphoton2, for file name
        num_frames: can come up with a better way to do this, currently just
         for knowing how many times to run loop2
        run_type: can be "image", "verify", "xylist", or "xylist_from_image".
         xylist runs DAOstarfinder and feeds list of sources to astrometry.net,
         image feeds the image directly to astrometry.net. image takes longest.
        dose: are we using dose maps?
        threshold: threshold for star detection used by dao, optional
        star_size: what size of star dao is looking for, will probably
         vary with the frame length used in the run
        output_directory: where to save the files
     """
    print(f"Running refiner for eclipse {eclipse}.")
    opt = {
        'eclipse': eclipse,
        'expt': expt,
        'num_frames': num_frames,
        'threshold': threshold,
        'star_size': star_size,
        'dose': dose,
        'run_type': run_type,
        'output_dir': output_directory
    }
    print(f"There are {num_frames} {expt} second depth frames for eclipse {eclipse}.")
    print(f"Making file name dictionary. File names will be saved to: {output_directory}")
    file_names = make_file_names(opt)
    # if you want to rerun the aspect correction to hopefully "refine" the solution, ends run
    # after verification
    if run_type == "xylist_from_image":
        # doesn't use pre-made xylists
        print("Running astrometry.net on xylists derived from frame images.")
        asp_df = run_xylist_on_image(eclipse, expt, num_frames, file_names, dose, **opt)
        print("Writing aspect solutions to csv.")
        asp_df.to_csv(file_names["asp_df"])
        return
    if run_type == "xylist":
        # fastest way to get aspect
        print("Running astrometry.net on xylists.")
        asp_df = run_xylist(eclipse, expt, num_frames, file_names, dose)
        print("Writing aspect solutions to csv.")
        asp_df.to_csv(file_names["asp_df"])
        return
    if run_type == "verify":
        # runs on images with a prior wcs to check it
        print("Running verification on previous aspect refinement.")
        verified_df = run_verification(eclipse, expt, num_frames, file_names, **opt)
        verified_df.to_csv(file_names["verified_df"])
        return
    if run_type == "image":
        # useful for testing the astrometry.net "in house" source extraction's performance
        print("Running astrometry.net directly on a frame images.")
        asp_df = run_image(eclipse, expt, num_frames, file_names)
        print("Writing aspect solutions to csv.")
        asp_df.to_csv(file_names["asp_df"])
        return
    print(f"Invalid run type: {run_type}.")
    return




