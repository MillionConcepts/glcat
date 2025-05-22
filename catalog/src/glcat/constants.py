"""
Named constants and types used throughout.
"""

from dataclasses import dataclass
from functools import total_ordering
from typing import Sequence

#: All the aperture sizes used during the GALEX mission.
ALL_APERTURES = [ 1.5, 2.3, 3.8, 6.0, 9.0, 12.8, 17.3, 30.]

#: Default set of aperture sizes to use for photometry.
DEFAULT_APERTURES = ALL_APERTURES

#: Default number of seconds to integrate over for each movie frame.
DEFAULT_DEPTH = 120


class BandMeta(type):
    """
    Metaclass for `Band`.
    This is a metaclass, rather than a __new__ hook on Band itself,
    because a __new__ hook cannot bypass __init__.
    """
    def __call__(cls, name, *args, **kwargs):
        if args or kwargs:
            if hasattr(cls, "NUV") and hasattr(cls, "FUV"):
                raise RuntimeError("creating a third band")
            return super(BandMeta, cls).__call__(name, **kwargs)

        match name.lower():
            case "nuv" | "near":
                return cls.NUV
            case "fuv" | "far":
                return cls.FUV
            case _:
                raise ValueError(f"unknown band {name!r}")


@dataclass(slots=True, frozen=True)
@total_ordering
class Band(metaclass=BandMeta):
    """
    A single band of the raw data.  There are exactly two instances of
    this class, Band.NUV and Band.FUV.  If 'xxx' is a recognized
    identifier for one of the two bands, Band("xxx") will return one
    of these instances.

    Properties are:
    name:       The official name of the band, either 'NUV' or 'FUV'.
    other:      The other band: Band.NUV.other == Band.FUV and vice versa.
    mag_scale:  Scale factor for conversion from counts per second in
                this band, to AB magnitude.
    flux_scale: Scale factor for conversion from counts per second in
                this band, to spectral flux density in cgs units
                (erg s⁻¹ cm⁻² Å⁻¹).
    """
    name: str
    other: "Band"
    mag_scale: float
    flux_scale: float

    def __init__(
        self,
        name: str,
        *,
        mag_scale: float,
        flux_scale: float,
    ):
        # have to use object.__setattr__ here to bypass the effect
        # of @dataclass(frozen=True)
        object.__setattr__(self, "name", name)
        object.__setattr__(self, "other", self)  # placeholder
        object.__setattr__(self, "mag_scale", mag_scale)
        object.__setattr__(self, "flux_scale", flux_scale)

    def __str__(self):
        return self.name

    def __repr__(self):
        return "Band." + self.name

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return False
        return self.name == other.name

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, type(self)):
            return NotImplemented
        # Reverse the usual ordering of strings, so that NUV < FUV.
        return self.name > other.name


Band.NUV = Band("NUV", mag_scale=20.08, flux_scale=2.06e-16)
Band.FUV = Band("FUV", mag_scale=18.82, flux_scale=1.4e-15)
# see comments in Band.__init__
object.__setattr__(Band.NUV, "other", Band.FUV)
object.__setattr__(Band.FUV, "other", Band.NUV)

def parse_bands(s: str) -> Sequence[Band]:
    """
    Parse a string as a Band value or comma-separated sequence
    of Band values; produce a sequence of Band values.

    Parsing is case insensitive. 'all' and 'both' are recognized
    as shorthand for 'nuv,fuv'.
    """
    val = set()

    for token in s.lower().split(','):
        token = token.strip()
        if token == "":
            pass
        elif token == "all" or token == "both":
            val.add(Band.NUV)
            val.add(Band.FUV)
        else:
            val.add(Band(token.strip()))
    return sorted(val)
