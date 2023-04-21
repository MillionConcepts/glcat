# psf fitting to measure changes in star size for aspect refinement

from astropy.io import fits
from astropy.nddata import NDData
from astropy.table import Table
from astropy.stats import sigma_clipped_stats
from photutils.psf import extract_stars, EPSFBuilder
import pandas as pd
import numpy as np
from compare_aspect.plots import centile_clip


def run_psf_compare(file_names):
    """ extract cutouts of stars from two different images using star
    locations id'd in a photometry file, filtered by area & flux to
    get relatively bright and small stars, fit a gaussian psf to each
    star. returns a table of stars with their gaussian fits for both
    images. """
    # use photom results from gphoton to get star locations to compare
    stars_tbl = get_good_stars(file_names)
    # og stars cutout
    og_stars = cutout_stars(stars_tbl, file_names['old_image_file'])
    # new stars cutout (same star locations)
    new_stars = cutout_stars(stars_tbl, file_names['new_image_file'])
    # pic of star cutouts
    make_psf_plots(og_stars, file_names)

    psf_comparison_tab = pd.DataFrame()
    for x in range(int(len(og_stars)*.1)):
        # currently only fitting psfs for 10% of stars in the sample...
        # subject to change, this saves time right now
        og_result_tab = fit_gaussian_prf(og_stars[x].data)
        new_result_tab = fit_gaussian_prf(new_stars[x].data)
        if len(og_result_tab) != 0:
            result1 = og_result_tab.to_pandas()
            result1["star"] = x
            result1["aspect_type"] = 'og'
            psf_comparison_tab = pd.concat([psf_comparison_tab,
                                            result1],
                                           ignore_index=True)
        if len(new_result_tab) != 0:
            result2 = new_result_tab.to_pandas()
            result2["star"] = x
            result2["aspect_type"] = 'new'
            psf_comparison_tab = pd.concat([psf_comparison_tab,
                                            result2],
                                           ignore_index=True)

    return psf_comparison_tab


def fit_gaussian_prf(image):
    """fit gaussian prf to a single star without a fixed sigma
    (initial guess is 2), using iteratively subtracted photometry"""
    from astropy.modeling.fitting import LevMarLSQFitter
    from astropy.stats import gaussian_sigma_to_fwhm
    from photutils.background import MADStdBackgroundRMS, MMMBackground
    from photutils.detection import IRAFStarFinder
    from photutils.psf import (DAOGroup, IntegratedGaussianPRF,
                               IterativelySubtractedPSFPhotometry)
    bkgrms = MADStdBackgroundRMS()
    std = bkgrms(image)
    iraffind = IRAFStarFinder(threshold=3.5 * std,
                              fwhm=4,
                              minsep_fwhm=0.01, roundhi=5.0, roundlo=-5.0,
                              sharplo=0.0, sharphi=2.0)
    daogroup = DAOGroup(2.0 * 4)
    mmm_bkg = MMMBackground()
    fitter = LevMarLSQFitter()
    # Circular Gaussian model integrated over pixels.
    gaussian_prf = IntegratedGaussianPRF(sigma=2)
    gaussian_prf.sigma.fixed = False
    photometry = IterativelySubtractedPSFPhotometry(finder=iraffind,
                                                    group_maker=daogroup,
                                                    bkg_estimator=mmm_bkg,
                                                    psf_model=gaussian_prf,
                                                    fitter=LevMarLSQFitter(),
                                                    niters=1,
                                                    fitshape=(11, 11))
    result_tab = photometry(image=image)

    return result_tab


def get_good_stars(file_names):
    """get xy list of stars in astropy table, currently using
    gphoton2 main for star detection"""
    photom_table = pd.read_csv(file_names['old_photom_file'])
    # these column names in the photom table only apply if the
    # tables were produced with the new image seg used in new
    # extraction branch of gphoton
    photom_table = photom_table[photom_table['area'] > 16]
    photom_table = photom_table[photom_table['area'] > 60]
    photom_table = photom_table[photom_table['segment_flux'] > .5]
    stars_tbl = Table()
    stars_tbl['x'] = photom_table['xcentroid']
    stars_tbl['y'] = photom_table['ycentroid']
    return stars_tbl


def cutout_stars(stars_tbl, imagefile):
    """cutout stars from full-depth image using star table"""
    # open full depth image to clip it (imagefile can be old or new aspect
    # image file path)
    image = fits.open(imagefile)
    data = image[1].data
    # background subtraction
    mean_val, median_val, std_val = sigma_clipped_stats(data, sigma=2.0)
    data -= median_val
    # change nans and inf to 0.0 for nddata purposes for psf fitting
    data = np.nan_to_num(data, copy=True, nan=0.0, posinf=0.0, neginf=0.0)
    # The extract_stars() function requires the input data
    # as an NDData object.
    nddata = NDData(data=data)
    # extract stars, can be done on multiple images
    stars = extract_stars(nddata, stars_tbl, size=21)
    return stars


def make_psf_plots(stars, file_names):
    """plot 25 of the star cutouts used for the psf fitting"""
    from numpy import random
    import matplotlib.pyplot as plt
    nrows = 5
    ncols = 5
    fig, ax = plt.subplots(nrows=nrows, ncols=ncols, figsize=(5, 5),
                           squeeze=True)
    ax = ax.ravel()
    for i in range(nrows * ncols):
        # currently only plotting the first 25 stars
        ax[i].imshow(centile_clip(stars[i]), origin='lower', cmap='viridis')
        ax[i].set_yticklabels([])
        ax[i].set_xticklabels([])
    plt.savefig(file_names['star_cutouts'])
    return


def fit_epsf(file_names, imagefile):
    """fit epsf to cutouts of individual stars, NOT USED RN
    in pipeline (doesn't return sigma fit in table)"""
    # code modified from photutils building an effective psf

    # get xy list of stars in astropy table
    stars_tbl = get_good_stars(file_names)
    # open full depth image to clip it (imagefile can be old or new aspect
    # image file path)
    image = fits.open(imagefile)
    data = image[1].data
    # background subtraction
    mean_val, median_val, std_val = sigma_clipped_stats(data, sigma=2.0)
    data -= median_val
    # change nans and inf to 0.0 for nddata purposes for psf fitting
    data = np.nan_to_num(data, copy=True, nan=0.0, posinf=0.0, neginf=0.0)
    # The extract_stars() function requires the input data as an NDData
    # object.
    nddata = NDData(data=data)
    # extract stars, can be done on multiple images
    stars = extract_stars(nddata, stars_tbl, size=21)
    # fit effective point spread function
    epsf_builder = EPSFBuilder(oversampling=4, maxiters=3,
                               progress_bar=False)
    epsf, fitted_stars = epsf_builder(stars)

    return epsf