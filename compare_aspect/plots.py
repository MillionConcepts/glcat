
import matplotlib.pyplot as plt
import pandas as pd


def make_plots(file_names):
    """main plot-making function for comparing aspect soln results"""

    # full depth image plots
    # side by side image comparison
    image_plot(file_names)
    # aspect (ra and dec) for both, plotted on top of each other

    # roll, plotted on top of each other

    # fwhm data plotted next to each other

    return


def image_plot(file_names):
    from astropy.io import fits
    fig, axs = plt.subplots(2)
    old_image = fits.open(file_names['old_image_file'])
    new_image = fits.open(file_names['new_image_file'])
    axs[0].imshow(centile_clip(old_image[1].data, (0, 99)),
                  interpolation='none', cmap='Greys_r', origin='lower')
    axs[1].imshow(centile_clip(new_image[1].data, (0, 99)),
                  interpolation='none', cmap='Greys_r', origin='lower')
    plt.savefig(file_names['image_comparison'])
    return


def ra_dec_plot(old_aspect, new_aspect, file_names):
    fig, axs = plt.subplots(2)
    axs[0].scatter(old_aspect["ra_center"], old_aspect["dec_center"], c="red",
                      s=0.7)
    axs[0].set_title("new aspect")
    axs[1].scatter(new_aspect["ra"], new_aspect["dec"], c="blue",
                      s=0.7)
    axs[1].set_title("old aspect")
    plt.savefig(file_names['ra_dec_plot'])
    return


def roll_plot(old_aspect, new_aspect, file_names):
    fig, axs = plt.subplots(2)
    axs[0].scatter(old_aspect.index, old_aspect["roll"], c="red",
                      s=0.7)
    axs[0].set_title(file_names['old_asp_title'])
    axs[1].scatter(new_aspect.index, new_aspect["orientation"]-360, c="blue",
                      s=0.7)
    axs[1].set_title(file_names['new_asp_title'])
    plt.savefig(file_names['roll_plot'])
    return


def centile_clip(image, centiles=(0, 90)):
    """
    simple clipping function that clips values above and below a given
    percentile range
    """
    import numpy as np

    finite = np.ma.masked_invalid(image)
    bounds = np.percentile(finite[~finite.mask].data, centiles)
    result = np.ma.clip(finite, *bounds)

    if isinstance(image, np.ma.MaskedArray):
        return result

    return result.data