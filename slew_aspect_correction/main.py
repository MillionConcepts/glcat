""" main pipeline for processing 1s of a slew frame """
from backplanes import make_backplanes
from astropy.io import fits
import numpy as np
import pandas as pd



def refine_slew_frame(frame_series):
    """ refine a slew frame, includes production of 10 .1 s backplanes """
    return


def slew_frame_pipeline(eclipse, frame_info):
    """ pipeline for a 1s slew frame. calls for 10 (.1 s) backplane images of the frame,
    uses a filtered 1s image to determine direction of movement via edge detection,
     stacks .1s images, and does source extraction. """
    print("filtering image")
    filter_image(frame_info)
    #TODO: run ASTRIDE using new file names
    print("running backplanes to make 1/10s frames")
    make_backplanes(
        eclipse=eclipse,
        band="NUV",
        depth=.1,
        leg=0,
        threads=4,
        burst=True,
        #TODO: change to use file management / new filenames
        local="/home/bekah/gphoton_working/test_data",
        kind="dose",
        radius=600,
        write={'array': True, 'xylist': False},
        inline=True,
        threshold=.45,
        star_size=2,
        snippet=frame_info['snippet']
    )
    streak_and_stack(frame_info)
    print("getting xylist for stacked slew frame")
    h, w = getXYlist_stackedFrames(frame_info)
    frame_info['h'] = h
    frame_info['w'] = w
    print("running astrometry on xylist for slew frame")
    astrometry_xylist_slew(frame_info)
    print("slew frame completed")
    return


def filter_image(frame_info):
    """use gaussian filter to smooth image for better processing"""
    from astropy.io import fits
    from skimage import filters
    dose_ais = fits.open(frame_info['1s_dose'])
    smooth = filters.gaussian(dose_ais[0].data, sigma=2)
    hdu = fits.PrimaryHDU(smooth)
    hdul = fits.HDUList([hdu])
    hdul.writeto(frame_info['smoothed_1s_dose'])
    return


def streak_and_stack(frame_info):
    """looks at ASTRIDE results and uses streak length / direction to
    stack .1s frames into a 1s image """
    #TODO: edit ones_frame var and layers
    expt = 2
    streaks = pd.read_csv(frame_info['streaks'])
    # calculate offsets
    streaks['x_diff'] = streaks['x_min'] - streaks['x_max']
    streaks['y_diff'] = streaks['y_min'] - streaks['y_max']
    # max offset
    ymax = max(abs(streaks['y_diff']))
    xmax = max(abs(streaks['x_diff']))
    ymean = np.mean(abs(streaks['y_diff']))
    xmean = np.mean(abs(streaks['x_diff']))
    slope = ymax / xmax
    mean = np.mean(streaks['slope'])
    print(f'calculated slope is {slope}, mean slope from streaks is {mean}.')
    xrate = xmax / expt
    yrate = ymax / expt  # used to be ymax / xmax
    # Number of pixels to offset each image.
    x_offset, y_offset = int(xrate * .1) - 5, int(yrate * .1) - 5
    frames = int((ones_frame) / .1)
    hdul = fits.open(f"/home/bekah/gphoton_working/test_data/e10982/e10982-nd-t00.1-b01-f{frames}-g_dose.fits.gz")
    img1 = hdul[0].data
    new_shape = ((layers - 1) * y_offset + img1.shape[0],
                 (layers - 1) * x_offset + img1.shape[1])
    stacked = np.zeros(new_shape)  # , dtype=np.float)
    stacked2 = np.zeros(new_shape)  # , dtype=np.float)
    # adding image layers together
    for layer in range(layers):
        print("adding image layers together")
        frame = frames + layer
        frame2 = frames + layer + layers
        dose_ais = fits.open(
            f"/home/bekah/gphoton_working/test_data/e10982/e10982-nd-t00.1-b01-f{frame}-g_dose.fits.gz")
        img1 = dose_ais[0].data  # *layer # for tagging pixels as being from a layer
        layer_op = (layers - 1) - layer
        stacked[layer_op * y_offset:layer_op * y_offset + img1.shape[0],
        layer_op * x_offset:layer_op * x_offset + img1.shape[1]] += img1
    # saving image to fits
    hdu = fits.PrimaryHDU(stacked)
    hdul = fits.HDUList([hdu])
    hdul.writeto(frame_info['stacked'], overwrite=True)

    return


def getXYlist_stackedFrames(filenames):
    """Use DAOStarFinder to extract point sources from stacked frames and
    produce a fits XYlist. """
    from photutils import DAOStarFinder
    # opening image for DAO starfinder
    s = fits.open(filenames['stacked'])
    sim = s[0].data
    # trying dao star finder / xylist call to astrometry.net instead of
    # using an image run
    daofind = DAOStarFinder(fwhm=4, threshold=1, sharplo=0.00)
    star_list = daofind(sim)
    tbl = star_list.to_pandas()
    star_list = tbl.sort_values(by="flux", ascending=False)
    colx = fits.Column(name='X', format='E', array=tbl['xcentroid'])
    coly = fits.Column(name='Y', format='E', array=tbl['ycentroid'])
    hdu = fits.BinTableHDU.from_columns([colx, coly])
    tableName = "starlist_slew.fits"
    hdu.writeto(tableName, overwrite=True)
    print('getting image width and height')
    h, w = sim.shape
    print(f'height is {h} and width is {w}')
    return h, w


def astrometry_xylist_slew(frame_info):
    import subprocess
    cmd = f"solve-field --overwrite -w {frame_info['w']} -e {frame_info['h']} " \
          f"--scale-units arcsecperpix --scale-low 1.0 --scale-high 1.5 " \
          f"--radius 5 {frame_info['stacked']}"
    subprocess.call(cmd, shell=True)
    print("Astrometry.net subprocess called.")
    return