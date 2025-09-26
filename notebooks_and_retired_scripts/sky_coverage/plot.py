EXPOSURE_LBLS = {
    0: "0",
    log10(60):    "1 m",
    log10(300):   "5 m",
    log10(900):   "15 m",
    log10(2700):  "45 m",
    log10(7980):  "2.2 h",
    log10(24000): "6.6 h",
    log10(72000): "20 h",
}

EXPOSURE_TICKS = sorted(EXPOSURE_LBLS.keys())
EXPOSURE_FMT = FuncFormatter(lambda x, _: EXPOSURE_LBLS[x])


def plot_exposure_from_stdin():
    with sys.stdin as fp:
        ra, dec, dose = np.genfromtxt(fp, skip_header=1, delimiter=",").T

    scaled_dose = np.zeros_like(dose)
    nonzero_dose = dose > 0
    scaled_dose[nonzero_dose] = np.log10(dose[nonzero_dose])

    fig, ax = plt.subplots(
        figsize=(20, 10),
        layout=ConstrainedLayoutEngine(
            h_pad=0, w_pad=0, hspace=0, wspace=0
        ),
        subplot_kw={
            "xlim": (-180, 180),
            "ylim": (-90, 90),
            "projection": ccrs.EqualEarth(globe=CEL_SPHERE)
        },
    )
    try:
        ctr = ax.hexbin(ra, dec, scaled_dose,
                        gridsize=3600,
                        transform=SKY_PLATE, cmap='Greys')
        fig.colorbar(
            ctr,
            location='top',
            fraction=0.025,
            pad=0.0125,
            ticks=EXPOSURE_TICKS,
            format=EXPOSURE_FMT,
            label="Exposure time",
        )
        fig.savefig("sky-coverage-full.png")
    finally:
        plt.close(fig)
