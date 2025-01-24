"""
Named constants and types used throughout.
"""

import enum

#: All the aperture sizes used during the GALEX mission
ALL_APERTURES = [ 1.5, 2.3, 3.8, 6.0, 9.0, 12.8, 17.3, 30., 60., 90., ]

#: Default set of aperture sizes to use for photometry.
DEFAULT_APERTURES = ALL_APERTURES

#: Default number of seconds to integrate over for each movie frame.
DEFAULT_DEPTH = 120


class Band(enum.Flag):
    """Indicates which band(s) of the raw data to process in a particular
       operation."""
    NUV = 0x01
    FUV = 0x02
    ALL = 0x03

    def __iter__(self) -> 'Iterator[Band]':
        """Iterate over all bits set in self."""
        if self & Band.NUV:
            yield Band.NUV
        if self & Band.FUV:
            yield Band.FUV

    @property
    def suffix(self):
        """Returns the movie file suffix for self.  Raises ValueError
           if self does not describe a single band."""
        if self == Band.NUV:
            return "mon"
        elif self == Band.FUV:
            return "mof"
        else:
            raise ValueError("suffix is only defined for single bands")

    @classmethod
    def parse(cls, s: str) -> 'Band':
        """Parse a string as a Band value or comma-separated sequence
           of Band values which are or-ed together.  Parsing is case
           insensitive and accepts several convenient aliases."""

        val = cls(0)
        for token in s.lower().split(','):
            token = token.strip()
            if token in ("nuv", "near"):
                val |= cls.NUV
            elif token in ("fuv", "far"):
                val |= cls.FUV
            elif token in ("both", "all"):
                val |= cls.ALL
            else:
                raise ValueError(f"unrecognized frequency band {token!r}")
        return val
