import numpy as np
import warnings

def counts2mag(cps, band):
    """
    Converts GALEX counts per second to AB magnitudes.
        See: http://asd.gsfc.nasa.gov/archive/galex/FAQ/counts_background.html

    :param cps: The flux in counts per second.

    :type cps: float

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :returns: float -- The converted flux in AB magnitudes.
    """

    scale = 18.82 if band == 'FUV' else 20.08

    # This threw a warning if the countrate was negative which happens when
    #  the background is brighter than the source. Suppress.
    with np.errstate(invalid='ignore'):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            mag = -2.5 * np.log10(cps) + scale

    return mag