There are three catalog tables generated per visit. The band-specific tables contain one entry / row per source detected in the visit-level image corresponding to the band of that table. So the NUV table (or `ncat`) contains one entry / row per source detected by DAOphot (via gPhoton2) in the visit-level NUV image. It then contains summary statistis for band-specific photometry as well as summary statistics for alternate-band photometry, when that band is valid, using identical sky coordinates for extraction on the alternate-band image. Another way to say this is that it contains "forced photometry" on the alternate band image. So the `ncat` contains photometry for NUV as well as FUV on the same source positions determined from the NUV image. ANd the `fcat` contains photometry for FUV as well as NUV on the same soure positions determined from the FUV image. Each band-specific table also contains a set of columns that describe the crossmatch between the `ncat` and `fcat`. The current columns are as follows.

### Band-Specific Catalog Column Definitions

| COLUMN NAME | DEFINITION | NOTES |
|-------------|------------|-------|
| `ECLIPSE` | The GALEX eclipse number of the observation. |
| `LEG` | **NOT YET IMPLEMENTED**
| `RA` / `DEC` | Right Ascension and Declination in J2000 decimal degrees at the center of the photometric aperture. | This is the location of a source detection as determined by DAOphot as implemented in the gPhoton2 software. |
| `[NF]UV_[XY]CENTER` | The x and y pixel coordinate in the corresponding band's image at the center of the photometric aperture. | The images have accurate WCS information, so they can be aligned in sky coordinates, but they do not generally (never?) have the same center coordinate or dimensions and therefore are not pre-aligned in image coordinates |
| `[NF]UV_SUM_APER[0-6]` | The sum of response-weighted counts inside of the photometric aperture in the particular band. | The `APER[0-6]` is a mapping to a specific set of photometric aperture radii that were used in the original mission-produced photometic catalogs. The mapping and related information appears in a table below. |
| `[NF]UV_EDGE_APER[0-6]` | **[TODO: Redundant. Remove from table.]**  |
| `[NF]UV_MASK_APER[0-6]` | **[TODO: Redundant. Remove from table.]**  |
| `[NF]UV_CPS_APER[0-6]` | The counts per second measured inside of the specified aperture. | This is the value of `[NF]UV_SUM_APER[0-6]` divided by the value of `[NF]UV_EXPT`. Note that because `[NF]UV_EXPT` is the entire visit integration depth, it does not correct for spatial variation in exposure due to detector motion. This should be disambiguated by using `[NF]UV_EDGE_FLAG_APER[0-6]` to determine whether the source crossed into the detector edge region (where it might then be affected by uncorrected spatial variations in exposure). The GLCAT measurements of such sources should be at least viewed with skepticism and possibly not used at all. If measurements of sources that are so affected are desired, then the correct approach is to use gPhoton2 to generate time-resolved observations and then build a coadd of frames in which the source is actually on the detector. |
| `[NF]UV_CPS_ERR_APER[0-6]` | The 1-sigma counting error on `[NF]UV_CPS_APER[0-6]`. | This is the square root of `[NF]UV_SUM_APER[0-6]` divided by `[NF]UV_EXPT` |
| `[NF]UV_MAG_APER[0-6]` | The AB Magnitude equivalent of `[NF]UV_CPS_APER[0-6]` using the conversion fact that `mag = -2.5 * np.log10(cps) + scale` where `scale = 20.08` for NUV and `scale = 18.82` for NUV. | An aperture correction has _not_ been applied to this value. Canonical aperture correction magnitudes can be found in the aperture definition table. |
| `[NF]UV_MAG_ERR_APER[0-6]` | The 1-sigma counting error `[NF]UV_MAG_APER[0-6]` expressed in the same units. | Note that counting errors are not _in general_ symmetrical in magnitude units, because of the logarithmic conversion. However they are approximately symmetrical when small. Also, small differences in magnitudes are approximately equivalent to percent differences, so this number can also be reasonably intepreted as a percentage in most cases. **[TODO: We should probably report upper and lower errors anyway.]** |
| `[NF]UV_MASK_FLAG_APER[0-6]` | A binary flag that indicates whether any events within the aperture were covered by by the hotspot mask. The nominal / unmasked value is `0`. |  |
| `[NF]UV_EDGE_FLAG_APER[0-6]` | A binary flag that indicates whether any events within the aperture were covered by by the edge mask. **[TODO: What distance from detector center are we currently using for the ede mask?] The nominal / unmasked value is `0`. |  |
| `[NF]UV_EXPT` | The effective exposure time of the visit. | "Effective" exposure time is the total observation time scaled for dead time effects and minus any periods of total data dropout. It is not corrected for spatial variation in exposure due to detector motion; see the note on `[NF]UV_CPS_APER[0-6]` about this. |
| `[NF]UV_EXTENDED` |  |  |
| ECLIPSE | **[TODO: Redundant. Remove from table.]** |  |
| `[FN]2[NF]_BEST_MATCH_INDEX`* | The index in the alternate-band table of the best sky-positional match (in angular distance) from that table to this one. | This does not necessarily indicate that they are the same source. See the note below on assessing crossmatch quality. This value can be empty if this is not the closest source to any source in the alternate-band catalog.\*\* This is most likely near the detector edges or in crowded fields. |
| `[FN]2[NF]_BEST_MATCH_DISTANCE` | The angular distance to the best sky-positional match to this source, as specififed by `[FN]2[NF]_BEST_MATCH_INDEX`. | This value can be empty. |
| `[FN]2[NF]_GOOD_MATCH_COUNT` | The number of sources in the alternate-band catalog within **[XX]** arcseconds of the source position.  | If this number is `0` then the "best match" is not also a "good match." If this number is `>1` then the field is almost certainly crowded and the crossmatch is subject to source confusion. If number is exactly `=1` then this _might_ be a good astrometric crossmatch to the source. Photometric contamination by nearby sources is still entirely possible outside of the `GOOD_MATCH` radius is still entirely possible. |
| [NF]2[FN]_MATCH_INDEX | The index in the alternate-band table of the best sky-positional match (in anular distance) from this band table to the alternate band table. | None of these values will be empty, as long as there are _any_ sources in the alternate-band table. |
| [NF]2[FN]_NUV_MATCH_INDEX | **[TODO: Vestigial. Remove.]** |  |
| [NF]2[FN]_MATCH_DISTANCE | The angular distance to the best sky-positional match from this band table to the alternate band table. |  |
| [NF]2[FN]_GOOD_MATCH_FLAG | This flag is set to `1` if the match from this band catalog to the alternate-band catalog is within **[XX]** arcseconds. | The inverse operation / summation of this flag is what produces the `GOOD_MATCH_COUNT` value in the alternate-band catalog. |
| [NF]2[FN]_BEST_MATCH_FLAG | This flag is set to `1` if the match from this band to the alternate-band catalog is also the best (i.e. closest) matching source to the alternate-band catalog. | Each entry in this catalog with this flag set will its index listed in the `BEST_MATCH_INDEX` column of the alternate-band catalog. |

\* The naming syntax is `[alternate-band initial]2[this-band initial]` and indicates that the first band has been matched _to_ the second band to produce this value. What this means is that every source in the first catalog has been assigned an index that matches to its sinle nearest match in the second catalog. Any particular source in the second catalog _might_ be the nearest match to _zero, one, or many_ sources in the first catalog.
\*\* The "empty value" constant is `-1`.

### Aperture Definition Table

| APER NUMBER | RADIUS (as) | NUV CORRECTION (AB Mag)* | NUV CORRECTION (AB Mag)* |
|-------------|-------------|--------------------------|--------------------------|
| 0 | 1.5 | foo | foo |
| 1 | 2.3 |  | foo |
| 2 | 3.8 |  | foo |
| 3 | 6.0 |  | foo |
| 4 | 9.0 |  | foo |
| 5 | 12.8 |  | foo |
| 6 | 17.3 |  | foo |

\* The aperture correction is not exactly correct for the whole detector.