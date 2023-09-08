import numpy as np
import pandas as pd
from astropy.io import fits
from backplanes import make_backplanes
from util import get_aspect_from_wcsinfo


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
    slew_paths = make_slew_paths(frame_series, aspect_root)
    # run ASTRIDE on 1s frame that already exists after filtering
    filter_image(frame_series, slew_paths)
    streaks = run_astride(frame_series, slew_paths)
    # if there are no detectable streaks, there's no point to this
    if streaks is None:
        return None
    # need to make 10 .1 s backplanes of the frame
    #TODO: check for backplane production
    make_shorter_backplanes(frame_series, eclipse_info, aspect_root)
    # stack frames based on streaks
    stack_frames(streaks, frame_series, slew_paths)
    # get xylist and send to astrometry.net
    h, w = getXYlist_stackedFrames(slew_paths)
    astrometry_xylist_slew(frame_series, slew_paths, h, w)
    aspect = get_aspect_from_wcsinfo(frame_series['wcs_path'])
    return aspect


def make_slew_paths(frame_series, aspect_root):
    slew_paths = {}
    slew_paths["smooth_backplane"] =
    slew_paths["short_backplanes"] =
    slew_paths["streaks"] =
    slew_paths["stacked"] =
    slew_paths["xylist"] =

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
        snippet=(eclipse_info['time'], eclipse_info['time']+1)
    )
    return


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

    streaks = get_streaks(slew_paths)

    return streaks


def get_streaks(slew_paths):
    """ get streaks from ASTRIDE results and calculate parameters """
    #TODO: change expt
    expt = 2
    streaks = pd.read_csv(slew_paths['streaks'])
    # calculate offsets
    streaks['x_diff'] = streaks['x_min'] - streaks['x_max']
    streaks['y_diff'] = streaks['y_min'] - streaks['y_max']
    # max offset
    ymax = max(abs(streaks['y_diff']))
    xmax = max(abs(streaks['x_diff']))
    #ymean = np.mean(abs(streaks['y_diff']))
    #xmean = np.mean(abs(streaks['x_diff']))
    slope = ymax / xmax
    mean = np.mean(streaks['slope'])
    xrate = xmax / expt
    yrate = ymax / expt  # used to be ymax / xmax
    # Number of pixels to offset each image.
    x_offset, y_offset = int(xrate * .1) - 5, int(yrate * .1) - 5
    streaks = {'x_offset': x_offset,
               'y_offset': y_offset,
               'slope':slope,
               'ymax':ymax,
               'xmax':xmax}
    return streaks


def stack_frames(streaks, frame_series, slew_paths):
    """looks at ASTRIDE results and uses streak length / direction to
    stack .1s frames into a 1s image """
    layers = 2 #TODO: change
    frames = int((ones_frame) / .1)
    hdul = fits.open(frame_series['backplane_path'])
    img1 = hdul[0].data
    new_shape = ((layers - 1) * streaks['y_offset'] + img1.shape[0],
                 (layers - 1) * streaks['x_offset'] + img1.shape[1])

    stacked = np.zeros(new_shape)  # , dtype=np.float)
    stacked2 = np.zeros(new_shape)  # , dtype=np.float)
    # adding image layers together
    #TODO: edit for correct backplane names in loop for short ones
    for layer in range(layers):
        print("adding image layers together")
        frame = frames + layer
        frame2 = frames + layer + layers
        dose_ais = fits.open(
            f"/home/bekah/gphoton_working/test_data/e10982/e10982-nd-t00.1-b01-f{frame}-g_dose.fits.gz")
        img1 = dose_ais[0].data  # *layer # for tagging pixels as being from a layer
        layer_op = (layers - 1) - layer
        stacked[layer_op * streaks['y_offset']:layer_op * streaks['y_offset'] + img1.shape[0],
        layer_op * streaks['x_offset']:layer_op * streaks['x_offset'] + img1.shape[1]] += img1
    # saving image to fits
    hdu = fits.PrimaryHDU(stacked)
    hdul = fits.HDUList([hdu])
    hdul.writeto('stacked2.fits', overwrite=True)
    return


def getXYlist_stackedFrames(slew_paths):
    """Use DAOStarFinder to extract point sources from stacked frames and
    produce a fits XYlist. """
    from photutils import DAOStarFinder
    s = fits.open(slew_paths['stacked_image'])
    sim = s[0].data
    daofind = DAOStarFinder(fwhm=4, threshold=1, sharplo=0.00)
    star_list = daofind(sim)
    tbl = star_list.to_pandas()
    #TODO: is output list sorted by flux? check
    star_list = tbl.sort_values(by="flux", ascending=False)
    colx = fits.Column(name='X', format='E', array=star_list['xcentroid'])
    coly = fits.Column(name='Y', format='E', array=star_list['ycentroid'])
    hdu = fits.BinTableHDU.from_columns([colx, coly])
    hdu.writeto(slew_paths['xylist'], overwrite=True)
    print('getting image width and height')
    h, w = sim.shape
    print(f'height is {h} and width is {w}')
    return h, w


def astrometry_xylist_slew(frame_series, slew_paths, h, w):
    """ run astrometry on a slew frame """
    import subprocess
    cmd = f"solve-field --overwrite --dir { frame_series['aspect_output']} " \
          f"-w {w} -e {h} --scale-units arcsecperpix --scale-low 1.0 --scale-high 1.5 " \
          f"-3 {frame_series['ra']} -4 { frame_series['dec']} " \
          f"--radius 5 {slew_paths['xylist']}"
    subprocess.call(cmd, shell=True)
    print("Astrometry.net subprocess called.")
    return
