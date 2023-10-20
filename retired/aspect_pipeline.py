# for combining aspect correction pipeline for slew frames and normal
# aspect corrected frames per eclipse
# in the ideal pipeline, the scst files are a giant parquet table and
# the backplanes are pre-produced
import sys
from slew_aspect_correction.main import slew_frame_pipeline
from retired.filenames import get_slew_frame_files, get_normal_frame_files, get_main_files
from astropy.io import fits
from pyarrow import parquet
from retired.xylist_runs import get_stars
from aspect_correction.combine_wcs import execute_combine
from retired.main_compare_asp import write_aspect2
from retired.xylist_runs import run_astrometry_net
from tools.retrieve import get_eclipse_data
from backplanes import make_backplanes, dosemaps_just_for_timestamps
sys.path.append('/home/bekah/gPhoton2')
from gPhoton.pipeline import execute_pipeline
from gPhoton.reference import PipeContext


def aspect_pipe(eclipse, band):
    """ handles solving aspect for all eclipses. """
    # how many legs does this eclipse have
    #TODO: load just by eclipse
    meta = parquet.read_table('/home/bekah/gPhoton2/gPhoton/aspect/metadata.parquet',
                              columns=['eclipse', 'legs']).to_pandas()

    # make file names using known frames
    file_names = get_main_files(eclipse, band, normal_frames, slew_frames)
    # are the files pre-produced on AWS S3 bucket backplane_test?
    try:
        #TODO: check if querying / get command for S3 costs money if no file there,
        #otherwise we don't want to hit it with every single eclipse
        get_eclipse_data(eclipse)
    except:
        print("Eclipse not in S3 bucket. Running gPhoton and producing backplanes locally.")
        # gotta make our own eclipse files, including backplanes
        try:
            execute_pipeline(
                eclipse,
                band,
                depth=120,
                # integer; None to deactivate (default None)
                threads=4,
                # where to both write output data and look for input data
                local_root="/home/bekah/gPhoton2/test_data",
                # auxiliary remote location for input data
                # remote_root="/mnt/s3",
                recreate=True,
                # list of floats; relevant only to lightcurve / photometry portion
                aperture_sizes=[12.8],
                # actually write image/movie products? otherwise hold in memory but
                # discard (possibly after performing photometry).
                write={"movie": False, "image": False},
                coregister_lightcurves=False,
                # photonpipe, moviemaker, None (default None)
                stop_after='photonpipe',
                photometry_only=False,
                # None, "gzip", "rice"
                compression="rice",
                # use array sparsification on movie frames?
                lil=True,
                # write movie frames as separate files
                burst=False,
                extended_photonlist=True,
                # aspect file, don't need to set unless need to use alt
                # file, 'aspect2.parquet'
                aspect="aspect2"
            )
        except:
            print(f"something didn't work :( for {eclipse} ")
            return
        # if photonlist is produced, try to make backplanes
        try:
            legs = metadata[metadata['eclipse'] == eclipse]['legs'].item()
            legs = [0] if legs == 0 else tuple(range(legs))
            for leg in legs:
                make_backplanes(
                    eclipse=eclipse,
                    band=band,
                    depth=1,
                    leg=leg,
                    threads=4,
                    burst=True,
                    local="/home/ubuntu/gPhoton2/test_data",
                    kind="dose",
                    radius=600,
                    write={'array': True, 'xylist': False},
                    inline=True,
                    threshold=.45,
                    star_size=2,
                    snippet=None
                )
        except:
            print(f"something didn't work :( for {eclipse} ")
            return
    print("Checking for expected backplane inventory.")
        frame_catalog = check_backplanes(eclipse, file_names)
    # frame_catalog says which backplanes may need to be produced, could add another point here to
    # produce those additional backplanes

    # get tables of normal frames and slew frames from big scst file
    normal_frames, slew_frames = get_frames(eclipse, band)

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


def check_backplanes(eclipse, file_names, ctx):
    """check for backplanes that should be in the eclipse folder based on
    extended aspect file"""
    import os
    ctx = PipeContext(
        eclipse,
        band,
        depth,
        "gzip",
        local,
        leg=leg,
        threads=threads,
        burst=burst,
        write=write,
        start_time=1000,
        stop_after=stop_after,
        snippet=snippet
    )
    time_stamps = dosemaps_just_for_timestamps(ctx, 400)
    # are all filenames in eclipse folder
    #TODO: put file name in file name function and change here
    files = os.listdir(f'home/bekah/gPhoton2/test_data/{eclipse}')

    # how to compare files? make a list of expected files and merge with actual files?
    # this approach is probably quicker than iterating through all of them
    frame_catalog =

    return frame_catalog


def run_astrometry_net(
        xylist_path,
        output_path,
        image_width,
        image_height,
        ra,
        dec,
        crpix_x,
        crpix_y):
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

    print(f"Running astrometry.net on frame.")
    cmd = f"solve-field --overwrite --no-plots --dir {output_path} -w {image_width} -e {image_height} " \
          f" --scale-units arcsecperpix --scale-low 1.0 --scale-high 1.6" \
          f" -N none -U none --temp-axy -B none -M none -R none " \
          f" -3 {ra} -4 {dec} --crpix-x {crpix_x} --crpix-y {crpix_y} --radius 5 {xylist_path}"
    # --no-tweak -3 {ra} -4 {dec} --scale-units arcsecperpix --scale-low 1.48 --scale-high 1.52
    # f" --verify '/home/bekah/glcat/astrometry/e{eclipse}/e{eclipse}-nd-{expt}s-0-f0000-rice.fits' --verify-ext 1 " \
    # -L 1.2 -H 1.8 -u app
    # RADIUS USED TO BE 3

    subprocess.call(cmd, shell=True)

    return None


