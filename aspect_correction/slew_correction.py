import numpy as np
import pandas as pd
from astropy.io import fits



def filter_image(files):
    """use gaussian filter to smooth image for better processing"""
    from astropy.io import fits
    from skimage import filters
    dose_ais = fits.open(f"/home/bekah/gphoton_working/test_data/e10982/e10982-nd-t0002-b01-f00{f}-g_dose.fits.gz")
    smooth = filters.gaussian(dose_ais[0].data, sigma=2)
    hdu = fits.PrimaryHDU(smooth)
    hdul = fits.HDUList([hdu])
    hdul.writeto('smoothed.fits')
    return


def streak_and_stack():
    """looks at ASTRIDE results and uses streak length / direction to
    stack .1s frames into a 1s image """
    expt = 2

    streaks = pd.read_csv(f"/home/bekah/glcat/smoothf{f}/streaks.csv")
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
    hdul.writeto('stacked2.fits', overwrite=True)
    return


def getXYlist_stackedFrames():
    """Use DAOStarFinder to extract point sources from stacked frames and
    produce a fits XYlist. """
    from photutils import DAOStarFinder
    # opening image for DAO starfinder
    s = fits.open('/home/bekah/glcat/stacked2.fits')
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


def astrometry_xylist_slew(h, w):
    import subprocess
    cmd = f"solve-field --overwrite -w {w} -e {h} " \
          f"--scale-units arcsecperpix --scale-low 1.0 --scale-high 1.5 " \
          f"--radius 5 '/home/bekah/glcat/starlist_slew.fits'"
    subprocess.call(cmd, shell=True)
    print("Astrometry.net subprocess called.")
    return
