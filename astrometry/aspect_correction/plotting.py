# for plotting different aspect solutions for different eclipses
from matplotlib import pyplot as plt
import pandas as pd


def main_plot(eclipse):
    from pyarrow import parquet

    # loading aspect table
    parq = parquet.read_table('/home/bekah/gphoton_working/gPhoton/aspect/aspect.parquet')
    aspect = parq.to_pandas()

    asprta = aspect[aspect["eclipse"] == eclipse].reset_index()

    eclipse_pad = str(eclipse).zfill(5)

    astrometry_xy_tycho_1s = pd.read_csv(f"/home/bekah/glcat/astrometry/aspect_correction/astrometry_xy_tycho_1s_{eclipse_pad}")
    astrometry_xy_tycho_5s = pd.read_csv(f"/home/bekah/glcat/astrometry/aspect_correction/astrometry_xy_tycho_5s_{eclipse_pad}")
    astrometry_xy_gaia_5s = pd.read_csv(f"/home/bekah/glcat/astrometry/aspect_correction/astrometry_xy_tycho_5s_{eclipse_pad}")

    image_file = f"/home/bekah/gphoton_working/test_data/e{eclipse_pad}/e{eclipse_pad}-nd-tfull-b00-image-r.fits"

    plot_ra_dec(eclipse, asprta, astrometry_xy_tycho_1s, astrometry_xy_tycho_5s, image_file)

    plot_roll(eclipse, asprta, astrometry_xy_tyconcho_1s, astrometry_xy_tycho_5s, image_file)

    image_plot(image_file)

    return


def image_plot(image_file):
    from astropy.io import fits
    # plot clipped image of pipeline aspect
    img = fits.open(image_file)
    plt.imshow(centile_clip(img[1].data, (0, 99)), cmap='Greys_r', origin='lower')
    #plt.savefig(f"{image_file}.jpg")

    return


def plot_roll(eclipse, asprta, astrometry_xy_tycho_1s, astrometry_xy_gaia_1s,
              image_file):
    from pyarrow import parquet

    fig, axs = plt.subplots(2, 2)

    axs[0, 0].scatter(asprta.index, asprta["roll"], s=0.7)
    axs[0, 0].set_title("OG Pipeline")

    axs[1, 0].scatter(astrometry_xy_tycho_1s.index, astrometry_xy_tycho_1s["orientation"]-360, s=0.7)
    axs[1, 0].set_title("astrometry_xy_tycho_1s")
    #axs[1, 0].sharex(axs[0, 0])

    axs[0, 1].scatter(astrometry_xy_gaia_1s.index, astrometry_xy_gaia_1s["orientation"]-360, s=0.7)
    axs[0, 1].set_title("astrometry_xy_tycho_5s")

    from astropy.io import fits
    # plot clipped image of pipeline aspect
    img = fits.open(image_file)
    axs[1, 1].imshow(centile_clip(img[1].data, (0, 99)), cmap='Greys_r', origin='lower')

    #axs[1, 1].scatter(astrometry_xy_gaia_5s.index, astrometry_xy_gaia_5s["orientation"]+360, s=0.7)
    axs[1, 1].set_title("full depth image")
    fig.tight_layout()

    plt.title(f"e{eclipse}")
    plt.savefig(f"roll_{eclipse}.jpg")

    return


def plot_ra_dec(eclipse, asprta, astrometry_xy_tycho_1s, astrometry_xy_gaia_1s,
                image_file):

    fig, axs = plt.subplots(2, 2)

    axs[0, 0].scatter(asprta["ra"], asprta["dec"], c="blue", s=0.7)
    axs[0, 0].set_title("OG Pipeline")

    axs[1, 0].scatter(astrometry_xy_tycho_1s["ra_center"], astrometry_xy_tycho_1s["dec_center"], c="red",
                      s=0.7)
    axs[1, 0].set_title("astrometry_xy_tycho_1s")
    #axs[1, 0].sharex(axs[0, 0])

    axs[0, 1].scatter(astrometry_xy_gaia_1s["ra_center"], astrometry_xy_gaia_1s["dec_center"], c="green",
                      s=0.7)
    axs[0, 1].set_title("astrometry_xy_tycho_5s")

    #axs[1, 1].scatter(astrometry_xy_gaia_5s["ra_center"], astrometry_xy_gaia_5s["dec_center"], c="orange",
    #                  s=0.7)
    from astropy.io import fits
    # plot clipped image of pipeline aspect
    img = fits.open(image_file)
    axs[1, 1].imshow(centile_clip(img[1].data, (0, 99)), cmap='Greys_r', origin='lower')
    axs[1, 1].set_title("full depth image")
    fig.tight_layout()

    plt.savefig(f"ra_dec_{eclipse}.jpg")

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