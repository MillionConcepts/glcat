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

def apcorrect1(radius, band):
    """
    Compute an apeture correction. First way. Uses the table data in Figure 4
        from Morissey, et al., 2007

    :param radius: The photometric radius, in degrees.

    :type radius: float

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :returns: float -- The aperture correction.
    """

    # [Future]: Handle arrays.

    if not band in ['NUV', 'FUV']:
        print("Invalid band.")
        return

    aper = np.array([1.5, 2.3, 3.8, 6.0, 9.0, 12.8, 17.3, 30., 60., 90.])/3600.

    if radius > aper[-1]:
        return 0.

    if band == 'FUV':
        dmag = [1.65, 0.96, 0.36, 0.15, 0.1, 0.09, 0.07, 0.06, 0.03, 0.01]
    else:
        dmag = [2.09, 1.33, 0.59, 0.23, 0.13, 0.09, 0.07, 0.04, -0.00, -0.01]
        if radius > aper[-2]:
            return 0.

    if radius < aper[0]:
        return dmag[0]

    ix = np.where((aper-radius) >= 0.)
    x = [aper[ix[0][0]-1], aper[ix[0][0]]]
    y = [dmag[ix[0][0]-1], dmag[ix[0][0]]]
    m, C = np.polyfit(x, y, 1)

    return m*radius+C