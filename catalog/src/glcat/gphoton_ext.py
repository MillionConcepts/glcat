"""
Material that should be merged into gphoton proper
"""

# Note: numpy-related type annotations in this file are not strictly
# correct, but are the best that can be done until changes happen in
# both core Python and numpy itself.  Further explanation here:
# https://numpy.org/doc/1.26/reference/typing.html#d-arrays

import numpy as np

from numpy.typing import ArrayLike
from glcat.constants import Band


def counts2mag(cps: ArrayLike, band: Band) -> np.ndarray:
    """
    Converts GALEX counts per second to AB magnitudes.
        See: http://asd.gsfc.nasa.gov/archive/galex/FAQ/counts_background.html

    :param cps: The flux in counts per second.

    :type cps: float

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :returns: float -- The converted flux in AB magnitudes.
    """

    # Negative count rates can happen due to the background being
    # contaminated by a nearby, brighter source.  Suppress warnings
    # from numpy; the result will be NaN, which is enough of an error
    # indicator.
    with np.errstate(invalid='ignore'):
        return -2.5 * np.log10(cps) + band.mag_scale


def counts2flux(cps: ArrayLike, band: Band) -> np.ndarray:
    """
    Converts GALEX counts per second to flux (erg sec^-1 cm^-2 A^-1).
        See: http://asd.gsfc.nasa.gov/archive/galex/FAQ/counts_background.html

    :param cps: The flux in counts per second.

    :type cps: float

    :param band: The band to use, either 'FUV' or 'NUV'.

    :type band: str

    :returns: float -- The converted flux in erg sec^-1 cm^-2 A^-1.
    """

    return np.asanyarray(cps) * band.flux_scale
