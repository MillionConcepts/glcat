"""
Loose methodology to stack slew frames together and get a refined
aspect solution:
1) identify slew frames by comparing time stamps of aspect parquet and
scst raw aspect solution (frames in scst file and not aspect parquet are
 considered slew frames if they have the correct voltage)
3) modify aspect parquet and run gphoton to make extended photonlist with
 slew frames
2) make 2s backplanes and .1 s backplanes made for the correct leg / eclipse
 / frames ID'd above
3) filter 2s backplane image to make it just blurry enough for ASTRIDE
to work
4) run ASTRIDE streak detection on the 2s backplanes to get streaks from
 image (ASTRIDE uses edge detection)
5) use streaks from ASTRIDE to get direction of movement (aka a linear best
 fit for streaks)
6) stack / combine 10 or 20 (still deciding) .1s frames to make a more
 "pointlike" image
7) use DAOstarfinder to get coordinates from stacked image and run
Astrometry.net on that xylist OR run Astrometry.net on the stacked image
8) backtrack to figure out which pixel in the stacked image is the "center"
 of the image and which timestamp it correlates with
9) edit aspect parquet with the astrometry.net solution
OPTIONAL:
10) run gphoton with edited aspect parquet
11) run gphoton with OG aspect soln from SCST file
12) check and compare FWHM of outputs from 10 and 11
"""

import pandas as pd
from astropy.io import fits
from pyarrow import parquet
import os
import numpy as np
import warnings
import sys
sys.path.insert(0, '/home/bekah/gphoton_working/gPhoton')


def nointernet_slew_pipeline(eclipse, band):
    """ for when there's no internet :( """
    # retrieve SCST, get info from eclipse header, and make df of slew frames
    scst, info = get_SCST(eclipse, band)
    # get slew frames
    slew_frames = get_slew_data(eclipse, scst, aspect_loc)
    # get slew frame chunks (frames no more than 1s apart, can assume from same leg)
    chunklist = get_chunks(slew_frames)
    # loop through slew stacking pipeline for multiple legs
    for i in range(len(chunklist)):
        print(f"Running slew stacking pipeline for slew frame chunk {i}.")
        slew_stacking_pipeline(chunklist[i], i, scst, info)
    print("Completed slew pipeline.")


def meta_slew_pipeline(eclipse, band, ):
    """for iterating over multiple chunks / legs of slew frames in an eclispe """

    aspect_loc = '/home/bekah/gphoton_working/gPhoton/aspect/aspect.parquet'

    # retrieve SCST, get info from eclipse header, and make df of slew frames
    scst, info = get_SCST(eclipse, band)
    # get slew frames
    slew_frames = get_slew_data(eclipse, scst, aspect_loc)
    # get slew frame chunks (frames no more than 1s apart, can assume from same leg)
    chunklist = get_chunks(slew_frames)
    # loop through slew stacking pipeline for multiple legs
    for i in range(len(chunklist)):
        print(f"Running slew stacking pipeline for slew frame chunk {i}.")
        slew_stacking_pipeline(chunklist[i], i, scst, info)
    print("Completed slew pipeline.")
    return


def slew_stacking_pipeline(frames, i, scst, info):
    """main pipeline for processing a chunk of slew frames from an eclipse.
     will have to be run multiple times for eclipses with multiple legs (ie AIS etc)
     that have slew frames between legs. """
    # generate file names
    filenames = make_file_names(i, info)
    # backplanes call

    #
    return


def make_file_names(i, info):

    files = {'aspect1': '/home/bekah/gphoton_working/gPhoton/aspect/aspect.parquet',
             'aspect2': ,
             'scst': ,
             'shortBackplanes': ,
             'longBackplanes': ,
             'smoothedLongBackplanes':
             'streaks': }
    return


def get_SCST(eclipse, band):
    """ Download SCST file from MAST """
    from gPhoton.io.mast import get_raw_paths, download_data
    paths = get_raw_paths(eclipse)
    scstpath = paths['scst']
    scst = download_data(
        eclipse, "scst", band, datadir=os.path.dirname(scstpath),
    )
    scst_pd = pd.DataFrame(scst[1].data)
    scst_pd = scst_pd.rename(columns={"pktime": "time"})

    # get eclipse info
    info = eclipse_info(scst)

    return scst_pd, info


def eclipse_info(scst):
    """ Get eclipse info: how many legs, length, flags, etc from SCST header """
    mtype = scst[0].header['MPSTYPE']
    legs = scst[0].header['MPSNPOS']
    ra_cent = scst[0].header['RA_CENT']
    dec_cent = scst[0].header['DEC_CENT']
    visit = scst[0].header['VISIT']
    einfo = {
             'mtype': mtype,
             'legs': legs,
             'ra_cent': ra_cent,
             'dec_cent': dec_cent,
             'visit': visit
             }
    return einfo


def get_slew_data(eclipse, scst_pd, aspect_loc):
    """
    Use SCST file to get missing slew data aka timestamps not available in
    the refined aspect soln that are in the scst file and are at HVNOM.
    """
    # loading aspect table
    parq = parquet.read_table(aspect_loc)
    aspect = parq.to_pandas()
    aspect = aspect[aspect["eclipse"] == eclipse]
    aspect = aspect.reset_index()

    # merge to get slew frames
    slew_frames = scst_pd.merge(aspect, how='left', on=['time'])
    slew_frames = slew_frames[slew_frames['hvnom_nuv'] == 1]

    return slew_frames


def get_chunks(slew_frames):
    """ returns chunklist, a dict of lists containing consecutive slew
    frame indexes from the scst table """
    indexlist = slew_frames.index.tolist()
    chunklist = {}
    counter = 0
    beginning = 0
    i = 0
    while len(indexlist) > i + 1:
        if indexlist[i] + 1 != indexlist[i + 1]:
            # then next number is not consecutive
            chunklist[counter] = indexlist[beginning:i + 1]
            beginning = i + 1
            counter = counter + 1
        i = i + 1
    return chunklist


def add_slew_to_aspect(slew_frames, files):
    """adds slewframes to aspect table for eclipse so that the slew
    frames can be added to photonlists"""
    # use acs solns for ra / dec / roll for slew frames in asp soln
    slew_frames['ra'] = slew_frames['ra'].fillna(slew_frames['ra_acs'])
    slew_frames['dec'] = slew_frames['dec'].fillna(slew_frames['dec_acs'])
    slew_frames['roll'] = slew_frames['roll'].fillna(slew_frames['roll_acs'])
    slew_frames = slew_frames.astype({
                          'ra': 'float64',
                          'dec': 'float64',
                          'roll': 'float64'})
    # save to parquet
    slew_frames.to_parquet(files["aspect2"], compression=None)
    print("modified aspect table to include slew frames.")
    return


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


# figure out any transformations required by center of image / pixel sitn
def correct_aspect():
    """getting ra, dec solution for the center of the first frame
    in the stacked series """

    return