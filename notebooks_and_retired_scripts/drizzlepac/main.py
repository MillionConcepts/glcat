from astropy.io import fits
from astropy.wcs import wcs
import drizzle
from matplotlib import pyplot as plt
import numpy as np
from astropy.io import ascii


def fix_wcs_per_frame():
    frames = ['0000', '0001', '0002', '0003']

    for frame in frames:
        hdulist = fits.open(
            f'/home/bekah/gphoton_working/test_data/e09869/files/e09869-nd-t0001-b00-f{frame}-g_dose.fits')
        header = hdulist[0].header
        f = fits.open(f'/home/bekah/glcat_tests/astrometry_temp/frame{frame[3]}_crop.wcs')
        w = WCS(f[0].header)
        hdulist[0].header.update(w.to_header())
        hdulist.writeto(
            f'/home/bekah/gphoton_working/test_data/e09869/files/e09869-nd-t0001-b00-f{frame}-  g_dose_WCS.fits',
            overwrite=True)
        fits.delval(f'/home/bekah/gphoton_working/test_data/e09869/files/e09869-nd-t0001-b00-f{frame}-g_dose_WCS.fits',
                    'CDELT1')
        fits.delval(f'/home/bekah/gphoton_working/test_data/e09869/files/e09869-nd-t0001-b00-f{frame}-g_dose_WCS.fits',
                    'CDELT2')
        fits.delval(f'/home/bekah/gphoton_working/test_data/e09869/files/e09869-nd-t0001-b00-f{frame}-g_dose_WCS.fits',
                    'PC1_1')
        fits.delval(f'/home/bekah/gphoton_working/test_data/e09869/files/e09869-nd-t0001-b00-f{frame}-g_dose_WCS.fits',
                    'PC1_2')
        fits.delval(f'/home/bekah/gphoton_working/test_data/e09869/files/e09869-nd-t0001-b00-f{frame}-g_dose_WCS.fits',
                    'PC2_1')
        fits.delval(f'/home/bekah/gphoton_working/test_data/e09869/files/e09869-nd-t0001-b00-f{frame}-g_dose_WCS.fits',
                    'PC2_2')
        hdulist.close()


def drizzle_demo_one(reference, outfile, infiles):
    """
    First demonstration of drizzle

    Parameters
    ==========
    reference
        A file containing the wcs of the output image

    outfile
        The name of the output image

    infiles
        The names of the input images to be combined
    """
    # Get the WCS for the output image
    hdulist = fits.open(reference)
    reference_wcs = wcs.WCS(hdulist[0].header)

    # Initialize the output with the WCS
    driz = drizzle.drizzle.Drizzle(outwcs=reference_wcs) #drizzle.

    # Combine the input images into on drizzle image
    for infile in infiles:
        infile = '/home/bekah/gphoton_working/test_data/e09869/files' + infile
        driz.add_fits_file(infile)

    # Write the drizzled image out
    driz.write(outfile)


def run_tweak_reg(infiles):
    from notebooks_and_retired_scripts.drizzlepac import tweakreg
    tweakreg.TweakReg(infiles,
                      use_custom_catalogs=True,
                      threshold=0.75,
                      conv_width=2,
                      # peakmax=None,
                      # imagefindcfg={'threshold' : 1000, 'conv_width' : 2},
                      # refimagefindcfg={'threshold' : 1000, 'conv_width' : 2},
                      # searchrad=3.0,
                      fitgeometry='general',
                      # exclusions='exclusions.txt',
                      computesig=True,
                      configobj=None,
                      interactive=False,
                      shiftfile=True,
                      outshifts='shift_default.txt',
                      updatehdr=False,
                      updatewcs=False)


def run_drizzle(infiles):
    # run astrodrizzle with a list of images
    from notebooks_and_retired_scripts.drizzlepac import astrodrizzle
    astrodrizzle.AstroDrizzle(
        input=('/home/bekah/gphoton_working/test_data/e09869/files/e09869-nd-t0001-b00-f0001-g_dose_WCS.fits',
        '/home/bekah/gphoton_working/test_data/e09869/files/e09869-nd-t0001-b00-f0002-g_dose_WCS.fits',
        '/home/bekah/gphoton_working/test_data/e09869/files/e09869-nd-t0001-b00-f0003-g_dose_WCS.fits'))

    astrodrizzle.AstroDrizzle(input=infile, output='final', wcskey='A', driz_sep_bits='64,32', final_wcs=True,
                              final_scale=1.5, final_rot=0)


def centile_clip(image, centiles=(1, 99)):
    """
    simple clipping function that clips values above and below a given
    percentile range
    """
    finite = np.ma.masked_invalid(image)
    bounds = np.percentile(finite[~finite.mask].data, centiles)
    result = np.ma.clip(finite, *bounds)

    if isinstance(image, np.ma.MaskedArray):
        return result

    return result.data


def plot_matches():
    match_tab = ascii.read('e09869-nd-t0001-b00-f0001-g_dose_WCS_catalog_fit.match')  # load match file in astropy table
    match_tab_chip2 = match_tab[match_tab['col15'] == 1]  # filter table for sources on chip 2 (on ext 1)
    x_cord, y_cord = match_tab_chip2['col11'], match_tab_chip2['col12']
    hdul = fits.open('/home/bekah/gphoton_working/test_data/e09869/files/e09869-nd-t0001-b00-f0001-g_dose_WCS.fits')
    plt.imshow(centile_clip(hdul[0].data, (0, 99.999)), interpolation=None)
    plt.scatter(x_cord, y_cord, s=30, edgecolor='r', facecolor='None', label='Matched Sources, Chip 2')
    # plt.ylim(0,2051)
    # plt.xlim(0,4096)
    # plt.legend(loc='best', fontsize=20)


def plot_stars_found():
    from matplotlib import pyplot as plt
    from astropy.table import Table

    # load xylist for 4 frames
    file0 = f'/home/bekah/gphoton_working/test_data/e09869/frame0.xyls'
    d0 = fits.open(file0)
    file1 = f'/home/bekah/gphoton_working/test_data/e09869/frame1.xyls'
    d1 = fits.open(file1)
    file2 = f'/home/bekah/gphoton_working/test_data/e09869/frame2.xyls'
    d2 = fits.open(file2)
    file3 = f'/home/bekah/gphoton_working/test_data/e09869/frame3.xyls'
    d3 = fits.open(file3)

    t0 = Table(d0[1].data)
    t1 = Table(d1[1].data)
    t2 = Table(d2[1].data)
    t3 = Table(d3[1].data)

    # plot x and y on frame 1 image
    hdul = fits.open('/home/bekah/gphoton_working/test_data/e09869/files/e09869-nd-t0001-b00-f0001-g_dose_WCS.fits')
    plt.imshow(centile_clip(hdul[0].data, (0, 99.999)), interpolation=None)
    plt.scatter(t0['X'], t0['Y'], c='red', alpha=.2)
    plt.scatter(t1['X'], t1['Y'], c='green', alpha=.2)
    plt.scatter(t2['X'], t2['Y'], c='orange', alpha=.2)
    plt.scatter(t3['X'], t3['Y'], c='blue', alpha=.2)
