# for combining aspect correction pipeline for slew frames and normal
# aspect corrected frames per eclipse
# in the ideal pipeline, the scst files are a giant parquet table and
# the backplanes are pre-produced

from slew_aspect_correction.main import slew_frame_pipeline
from filenames import get_slew_frame_files, get_normal_frame_files, get_main_files
from astropy.io import fits
from pyarrow import parquet
from aspect_correction.xylist_runs import get_stars
from aspect_correction.combine_wcs import execute_combine
from compare_aspect.main import write_aspect2
from aspect_correction.xylist_runs import run_astrometry_net


def aspect_pipe(eclipse, band):
    """ handles solving aspect for all eclipses. """
    # get tables of normal frames and slew frames from big scst file
    normal_frames, slew_frames = get_frames(eclipse, band)

    # make file names using known frames
    file_names = get_main_files(eclipse, band, normal_frames, slew_frames)
    # download backplanes from aws (?)
    #backplanes = get_backplanes()

    # loop through every aspect frame
    # this is where it could be easily parallelized
    print("Solving for aspect soln of normal frames.")
    for nframe in normal_frames['frame']:
        frame_info = get_normal_frame_files(nframe)
        #TODO: download frame from s3
        movie = fits.open(frame_info['1s_dose'])
        table_name = get_stars(movie[0].data, i, file_names)
        #TODO: change inputs in run_astrometry_net to work with frame_info
        run_astrometry_net(frame_info, table_name)
    print("Completed normal frames.")

    if len(slew_frames) != 0:
        print("Moving on to slew frames.")
        for sframe in slew_frames['frame']:
            frame_info = get_slew_frame_files(sframe)
            # TODO: download frame from s3, download extended photonlist (need
            #  for backplanes w/o gphoton run)
            slew_frame_pipeline(eclipse, frame_info)
        print("Completed slew frames.")

    print("Combining wcs for each frame.")
    #TODO: edit execute combine to work with new file names
    # (won't need expt & number of frames)
    wcs_df = execute_combine(eclipse,
                             expt,
                             num_frames,
                             "image_to_xylist",
                             dose=True,
                             output_directory="/home/bekah/glcat_tests/")

    aspect_file = write_aspect2(eclipse, file_names)
    print(f"Succesfully save aspect file to {aspect_file}, cleaning up wcs files.")

    # TODO: remove unnecessary files from eclipse folder
    #  (leave aspect soln, extended parquet, etc)
    #clean_up_wcs(file_names)

    return aspect_file


def get_frames(eclipse, band):
    """ look at aspect table and determine which frames are normal and slew frames """
    #TODO: the eclipse value can't be padded for filtering
    normal_frames = parquet.read_table('/home/bekah/gphoton_working/extended_scst_7_16.parquet',
                  filters=[('eclipse', '=', eclipse), ('slew', '=', False)]).to_pandas()

    slew_frames = parquet.read_table('/home/bekah/gphoton_working/extended_scst_7_16.parquet',
                  filters=[('eclipse', '=', eclipse), ('slew', '=', True)]).to_pandas()

    return normal_frames, slew_frames



