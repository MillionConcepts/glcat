import numpy as np
import pandas as pd
from astropy.io import fits
from backplanes import make_backplanes
from aspect_correction.util import get_aspect_from_wcsinfo
import os

def refine_slew_frame(
        frame_series,
        eclipse_info,
        aspect_root):
    """ refine a slew frame. steps:
    1) filter 1s backplane w gaussian
    2) identify streaks in filtered image
    3) create .1s backplanes from time snippet of frame
    4) stack .1s backplanes along streak
    5) use DAOstarfinder to get stars from stacked image
    6) astrometry.net on resulting xylist from 5)
    7) get aspect via wcs info"""
    # make paths for slew specific outputs
    slew_paths = make_slew_paths(frame_series)
    # run ASTRIDE on 1s frame that already exists after filtering
    filter_image(frame_series, slew_paths)
    streaks = run_astride(frame_series, slew_paths)
    # if there are no detectable streaks, there's no point to this
    if streaks is None:
        return None
    # need to make 10 .1 s backplanes of the frame
    # TODO: check for backplane production
    make_shorter_backplanes(frame_series, eclipse_info, aspect_root)
    if check_for_short_backplanes(slew_paths["short_backplanes"]):
        # get streak direction
        # direction = get_stack_direction(frame_series, slew_paths)
        # stack frames based on streaks
        stack_frames(streaks, frame_series, slew_paths)
        # get xylist and send to astrometry.net
        h, w = get_xylist_for_stacked(slew_paths, frame_series)
        astrometry_xylist_slew(frame_series, h, w)
        aspect = get_aspect_from_wcsinfo(frame_series['wcs_path'])
        return aspect
    print("Short backplane creation failed for this frame.")
    return


def make_slew_paths(frame_series):
    """ making paths for files that only exist for slew frames """
    slew_paths = {}
    slew_paths["smooth_backplane"] = frame_series["aspect_output"] + \
                                     f"smooth_{frame_series['time_stamp']}.fits"
    slew_paths["short_backplanes"] = [frame_series['backplane_path'] \
                                          .replace('f0001', 'f00001') for x in range(10)]
    slew_paths["streaks"] = frame_series["aspect_output"] + "streaks.txt"
    slew_paths["stacked"] = frame_series["aspect_output"] + \
                            f"stacked_{frame_series['time_stamp']}.fits"
    slew_paths["stacked_opposite"] = frame_series["aspect_output"] + \
                                     f"stacked_op_{frame_series['time_stamp']}.fits"
    return slew_paths


def make_shorter_backplanes(frame_series, eclipse_info, aspect_root):
    """ make backplanes of .1 s from 1s snippet that is a frame  """
    print("writing slew frame .1s backplanes")
    make_backplanes(
        eclipse=eclipse_info['eclipse'],
        band=eclipse_info['band'],
        depth=.1,
        leg=frame_series['leg'],
        threads=4,
        burst=True,
        local=aspect_root,
        kind="dose",
        radius=400,
        write={'array': True, 'xylist': False},
        inline=True,
        threshold=.75,
        star_size=2,
        snippet=(eclipse_info['time'], eclipse_info['time'] + 1)
    )
    return


def check_for_short_backplanes(short_backplanes):
    """ check that short backplanes were created
    before trying to stack """
    return all(list(map(os.path.isfile, short_backplanes)))


def filter_image(frame_series, slew_paths):
    """use gaussian filter to smooth image for better processing"""
    from astropy.io import fits
    from skimage import filters
    dose_ais = fits.open(frame_series["backplane_path"])
    smooth = filters.gaussian(dose_ais[0].data, sigma=2)
    hdu = fits.PrimaryHDU(smooth)
    hdul = fits.HDUList([hdu])
    hdul.writeto(slew_paths["smooth_backplane"])
    return


def run_astride(frame_series, slew_paths):
    """ run external program astride """
    # Import the ASTRiDE library.
    from astride import Streak
    # Read a fits image and create a Streak instance.
    streak = Streak(
        slew_paths["smooth_backplane"],
        contour_threshold=4,
        min_points=1800,
        radius_dev_cut=0.3,
        shape_cut=0.1)
    # Detect streaks.
    streak.detect()
    # Write outputs and plot figures.
    streak.write_outputs()
    streak.plot_figures()
    ### streak.streaks is
    # 1) an insane data structure that I don't understand why anyone would
    # voluntarily store their data that way!!!!!!!!!
    # 2) a list of dictionaries where each dictionary is a streak, and 'x'
    # 'y' within each dictionary each hold an array of x and y edge points,
    # respectively
    x_vals = []
    y_vals = []
    for s in streak.streaks:
        x_vals.append(s['x'])
        y_vals.append(s['y'])
    # these are for getting the outlines of streaks, which we don't really need
    # except for plotting
    # np.save('xvals.npy', np.array(x_vals, dtype=object), allow_pickle=True)
    # np.save('yvals.npy', np.array(y_vals, dtype=object), allow_pickle=True)
    streaks = get_streaks(slew_paths)
    return streaks


def get_streaks(slew_paths):
    """ get streaks from ASTRIDE results and calculate parameters """
    # TODO: change expt
    expt = 2
    streaks = pd.read_csv(slew_paths['streaks'], delim_whitespace=True)
    # calculate offsets
    streaks['x_diff'] = streaks['x_min'] - streaks['x_max']
    streaks['y_diff'] = streaks['y_min'] - streaks['y_max']
    # max offset
    ymax = max(abs(streaks['y_diff']))
    xmax = max(abs(streaks['x_diff']))
    slope = ymax / xmax
    # mean = np.mean(streaks['slope'])
    xrate = xmax / expt
    yrate = ymax / expt  # used to be ymax / xmax
    # Number of pixels to offset each image.
    x_offset, y_offset = int(xrate * .1) - 5, int(yrate * .1) - 5
    streaks = {'x_offset': x_offset,
               'y_offset': y_offset,
               'slope': slope,
               'ymax': ymax,
               'xmax': xmax}
    return streaks


def stack_frames(streaks, frame_series, slew_paths):
    """ stack frames in two directions along linear streak """
    # num frames
    num_frames = 10
    # dimensions of first frame
    first_frame_hdul = fits.open(frame_series['backplane_path'])
    first_frame = first_frame_hdul[0].data
    # new image size is the offset between frames * the number of frames being added
    # plus the og image size
    # num # frames is 10 usually (.1s increments for 1s)
    x_dim = ((num_frames - 1) * streaks['x_offset'] + first_frame.shape[0])
    y_dim = ((num_frames - 1) * streaks['y_offset'] + first_frame.shape[1])

    # empty image of stacked dimensions
    stacked = np.zeros((x_dim, y_dim))
    stacked_opposite = np.zeros((x_dim, y_dim))

    # iterate through adding 10 .1 s frames
    for frame in range(num_frames):
        # get frame image
        frame_hdul = fits.open(slew_paths['short_backplanes'][frame])
        frame_image = frame_hdul[0].data
        # countdown layers
        inverse_layer = (num_frames - 1) - frame
        # add frame to section of stacked image that corresponds with offset
        stacked[inverse_layer * streaks['y_offset']:inverse_layer * streaks['y_offset'] +
                                                    frame_image.shape[0],
        inverse_layer * streaks['x_offset']:inverse_layer * streaks['x_offset'] +
                                            frame_image.shape[1]] += frame_image
        # add frame to section of stacked image that corresponds with offset in opposite dir
        stacked_opposite[frame * streaks['y_offset']:frame * streaks['y_offset'] +
                                                     frame_image.shape[0],
        frame * streaks['x_offset']:frame * streaks['x_offset'] +
                                    frame_image.shape[1]] += frame_image
    print("Finished stacking.")
    # saving image to fits
    hdu = fits.PrimaryHDU(stacked)
    hdul = fits.HDUList([hdu])
    hdul.writeto(slew_paths["stacked"], overwrite=True)
    hdu2 = fits.PrimaryHDU(stacked_opposite)
    hdul2 = fits.HDUList([hdu2])
    hdul2.writeto(slew_paths["stacked_opposite"], overwrite=True)
    print("Wrote two stacked images.")

    return


def get_xylist_for_stacked(slew_paths, frame_series):
    """Use DAOStarFinder to extract point sources from stacked frames and
    produce a fits XYlist. """
    from photutils import DAOStarFinder

    s = fits.open(slew_paths['stacked_image'])
    sim = s[0].data
    daofind = DAOStarFinder(fwhm=4, threshold=1, sharplo=0.00)
    star_list = daofind(sim)

    s = fits.open(slew_paths['stacked_opposite'])
    sim = s[0].data
    daofind = DAOStarFinder(fwhm=4, threshold=1, sharplo=0.00)
    star_list_opposite = daofind(sim)

    # the logic here is that the correctly stacked image will
    # have less stars found by DAOStarFinder (less streaks)
    if len(star_list_opposite) > len(star_list):
        tbl = star_list.to_pandas()
    if len(star_list_opposite) <= len(star_list):
        tbl = star_list_opposite.to_pandas()

    star_list = tbl.sort_values(by="flux", ascending=False)
    colx = fits.Column(name='X', format='E', array=star_list['xcentroid'])
    coly = fits.Column(name='Y', format='E', array=star_list['ycentroid'])
    hdu = fits.BinTableHDU.from_columns([colx, coly])
    hdu.writeto(frame_series['xylist_path'], overwrite=True)
    print('getting image width and height')
    h, w = sim.shape
    print(f'height is {h} and width is {w}')
    return h, w


def astrometry_xylist_slew(frame_series, h, w):
    """ run astrometry on a slew frame """
    import subprocess
    cmd = f"solve-field --overwrite --dir {frame_series['aspect_output']} " \
          f"-w {w} -e {h} --scale-units arcsecperpix --scale-low 1.0 --scale-high 1.5 " \
          f"-3 {frame_series['ra']} -4 {frame_series['dec']} " \
          f"--radius 5 {frame_series['xylist_path']}"
    subprocess.call(cmd, shell=True)
    print("Astrometry.net subprocess called.")
    return


# **************** not used rn **********************************************************


def get_stack_direction(frame_series, slew_paths):
    """ get direction of slew, so we know which direction to stack short frames.
     uses og aspect, although you could use new if available """
    dir = frame_series['time_stamp'] + 1
    return dir
