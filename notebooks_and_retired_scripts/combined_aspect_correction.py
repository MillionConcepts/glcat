# for combining aspect correction pipeline for slew frames and normal
# aspect corrected frames per eclipse
# in the ideal pipeline, the scst files are a giant parquet table and
# the backplanes are pre-produced


# pipeline per eclipse


def main_aspect_pipe():
    """ handles solving aspect for all eclipses. """
    # get tables of normal frames and slew frames from big scst file
    normal_frames, slew_frames = get_frames()
    # download backplanes from aws (?)
    backplanes = get_backplanes()
    # loop through every aspect frame
    # this is where it could be easily parallelized
    print("Solving for aspect soln of normal frames.")
    for nframe in normal_frames:
        # want to get xylist in big run or get xylist indiv?
        # start with all xylists and then modify the run if necessary?
        frame_info = {}
        run_astrometry_net(frame_info)
        wcs = get_wcs(frame_info)
        # append wcs to df with frame info
    if len(slew_frames) != 0:
        print("Moving on to slew frames.")
        for sframe in slew_frames:
            frame_info = {}
            run_slew_frame(frame_info)


def run_astrometry_net(frame_info):
    import subprocess
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

    cmd = f"solve-field --overwrite --no-plots --dir {file_names['astrometry_temp']} -w {image_width} -e {image_height} " \
          f"--scale-units arcsecperpix --scale-low 1.48 --scale-high 1.52 " \
          f"-N none -U none --axy none -B none -M none -R none " \
          f" -3 {ra} -4 {dec} --crpix-x {crpix_x} --crpix-y {crpix_y} --radius 5 {file_names[xylist_type][i]}"

    subprocess.call(cmd, shell=True)
    return None


def get_wcs(frame_info):
    # if you keep the timestamp, it doesn't matter when they're added to the aspect table, it can be
    # sorted by time later
    return wcs_info