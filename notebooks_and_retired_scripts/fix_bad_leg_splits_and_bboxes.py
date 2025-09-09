"""Fix two classes of errors made by boresight_from_scst.py:

- Correct bounding box records for eclipses that cross the RA±180
  meridian.  Fully automated.

- Merge all the legs of eclipses that shouldn't have been split into
  multiple legs.  Manual: you must provide a list of eclipses to merge
  on the command line.
"""

import argparse
import os
import sys
import tempfile

from itertools import repeat
from math import pi as PI, inf as Inf
from pathlib import Path
from typing import Any, Iterable

import pyarrow.compute as pc
from pyarrow import Table, concat_tables, parquet

import shapely


SC = parquet.SortingColumn

EXPECTED_COLS_METADATA = frozenset((
    "eclipse",
    "plan_type",
    "plan_subtype",
    "visit",
    "plan_id",
    "planned_legs",
    "observed_legs",
    "has_aspect",
    "eclipse_start",
    "eclipse_duration",
    "ok_exposure_time",
    "ra0",
    "ra_min",
    "ra_max",
    "dec0",
    "dec_min",
    "dec_max",
    "nuv_det_on_time",
    "nuv_det_temp",
    "nuv_tdc_temp",
    "nuv_has_raw6",
    "fuv_det_on_time",
    "fuv_det_temp",
    "fuv_tdc_temp",
    "fuv_has_raw6",
))

EXPECTED_COLS_BORESIGHT = frozenset((
    "eclipse",
    "leg",
    "time",
    "duration",
    "ra0",
    "ra_min",
    "ra_max",
    "dec0",
    "dec_min",
    "dec_max",
    "planned_ra",
    "planned_dec",
    "full_exposure_area",
    "full_exposure_uncircularity",
    "partial_exposure_area_ratio",
))

EXPECTED_COLS_LEG_APERTURE = frozenset((
    "eclipse",
    "leg",
    "partial",
    "full_400",
    "full_370",
    "full_350",
))

EXPECTED_COLS_AFTER_JOINS = frozenset((
    "eclipse",
    "plan_type",
    "plan_subtype",
    "visit",
    "plan_id",
    "planned_legs",
    "observed_legs",
    "has_aspect",
    "eclipse_start",
    "eclipse_duration",
    "ok_exposure_time",
    "ra0_e",
    "ra_min_e",
    "ra_max_e",
    "dec0_e",
    "dec_min_e",
    "dec_max_e",
    "nuv_det_on_time",
    "nuv_det_temp",
    "nuv_tdc_temp",
    "nuv_has_raw6",
    "fuv_det_on_time",
    "fuv_det_temp",
    "fuv_tdc_temp",
    "fuv_has_raw6",
    "leg",
    "time",
    "duration",
    "ra0_l",
    "ra_min_l",
    "ra_max_l",
    "dec0_l",
    "dec_min_l",
    "dec_max_l",
    "planned_ra",
    "planned_dec",
    "full_exposure_area",
    "full_exposure_uncircularity",
    "partial_exposure_area_ratio",
    "partial",
    "full_400",
    "full_370",
    "full_350"
))


def check_columns(tbl: Table, expected: frozenset) -> None:
    cols = frozenset(tbl.column_names)
    if cols != expected:
        missing = ", ".join(sorted(expected - cols))
        extra = ", ".join(sorted(cols - expected))
        raise RuntimeError(
            f"column mismatch:\n"
            f" missing: {missing}\n"
            f"   extra: {extra}\n"
        )


def add_360_to_ra(coords):
    c = coords.copy()
    c[:,0] += 360.0
    return c


def recalculate_bounds(region):
    if isinstance(region, shapely.Polygon):
        ra_min, ra_max, dec_min, dec_max = region.bounds
        assert not (ra_min == -180.0 and ra_max == 180.0)
        ra0, dec0 = region.centroid.coords[0]
    else:
        assert isinstance(region, shapely.MultiPolygon)
        adjusted = shapely.union_all([
            (
                shapely.transform(r, add_360_to_ra)
                if r.bounds[0] == -180.0
                else r
            ) for r in region.geoms
        ])
        ra_min, ra_max, dec_min, dec_max = adjusted.bounds
        assert not (ra_min == -180.0 and ra_max == 180.0)
        ra0, dec0 = adjusted.centroid.coords[0]

    return ra0, ra_min, ra_max, dec0, dec_min, dec_max


def recalculate_metrics(partial, full):
    if full.is_empty:
        return (0, Inf, Inf)

    full_area = full.area
    if isinstance(full, shapely.Polygon):
        full_perimeter = full.exterior.length
    else:
        assert isinstance(full, shapely.MultiPolygon)
        full_perimeter = sum(p.exterior.length for p in full.geoms)

    # The standard isoperimetric quotient is 4πA⁄P², which ranges from
    # 0 to 1, where 1 is a perfect circle and 0 is a shape with zero
    # area but nonzero perimeter, e.g. a line segment.  Its reciprocal
    # ranges from 1 to positive infinity, and subtracting 1 puts a
    # perfect circle at 0.  We do this transformation because most
    # legs have an isoperimetric quotient very close to 1; on our data
    # set, the shifted reciprocal quotient ranges from 4 × 10⁻⁵ to
    # 0.21, which is usefully put on a log scale.)
    uncircularity = full_perimeter * full_perimeter / (4 * PI * full_area) - 1

    # For the "partial area ratio", what we want is the _exclusive_
    # partial area, that is, the area of the partial region minus the
    # area of the full region.
    partial_area_ratio = (partial.area - full_area) / full_area
    return (full_area, uncircularity, partial_area_ratio)


def merge_legs_and_recalculate_bounds(er):
    # these columns all come from the metadata table and their values
    # are all fine as is; we just need to squash them down to one row
    del er["eclipse"][1:]
    del er["plan_type"][1:]
    del er["plan_subtype"][1:]
    del er["visit"][1:]
    del er["plan_id"][1:]
    del er["planned_legs"][1:]
    del er["has_aspect"][1:]
    del er["eclipse_start"][1:]
    del er["eclipse_duration"][1:]
    del er["ok_exposure_time"][1:]
    del er["nuv_det_on_time"][1:]
    del er["nuv_det_temp"][1:]
    del er["nuv_tdc_temp"][1:]
    del er["nuv_has_raw6"][1:]
    del er["fuv_det_on_time"][1:]
    del er["fuv_det_temp"][1:]
    del er["fuv_tdc_temp"][1:]
    del er["fuv_has_raw6"][1:]

    # we now only have one leg
    er["observed_legs"][0] = 1
    del er["observed_legs"][1:]
    del er["leg"][1:]
    del er["planned_ra"][1:]
    del er["planned_dec"][1:]

    # the new start time is the same as the start time of the old leg 1;
    # because the leg now includes a bunch of intermediate flagged points,
    # the new duration is (old leg N's start time - old leg 1's start time)
    # + old leg N's duration.
    er["duration"][0] = (er["time"][-1] - er["time"][0]) + er["duration"][-1]
    del er["time"][1:]
    del er["duration"][1:]

    # reconstitute and merge the shapes
    partial = shapely.union_all([
        shapely.from_wkb(blob) for blob in er["partial"]
    ])
    er["partial"][0] = partial.wkb
    del er["partial"][1:]

    full_400 = shapely.intersection_all([
        shapely.from_wkb(blob) for blob in er["full_400"]
    ]).buffer(0)
    er["full_400"][0] = full_400.wkb
    del er["full_400"][1:]

    full_370 = shapely.intersection_all([
        shapely.from_wkb(blob) for blob in er["full_370"]
    ]).buffer(0)
    er["full_370"][0] = full_370.wkb
    del er["full_370"][1:]

    full_350 = shapely.intersection_all([
        shapely.from_wkb(blob) for blob in er["full_350"]
    ]).buffer(0)
    er["full_350"][0] = full_350.wkb
    del er["full_350"][1:]

    ra0, ra_min, ra_max, dec0, dec_min, dec_max = \
        recalculate_bounds(full_400 if partial.is_empty else partial)

    er["ra0_e"][0] = ra0
    er["ra0_l"][0] = ra0
    del er["ra0_e"][1:]
    del er["ra0_l"][1:]

    er["ra_min_e"][0] = ra_min
    er["ra_min_l"][0] = ra_min
    del er["ra_min_l"][1:]
    del er["ra_min_e"][1:]

    er["ra_max_e"][0] = ra_max
    er["ra_max_l"][0] = ra_max
    del er["ra_max_l"][1:]
    del er["ra_max_e"][1:]

    er["dec0_e"][0] = dec0
    er["dec0_l"][0] = dec0
    del er["dec0_l"][1:]
    del er["dec0_e"][1:]

    er["dec_min_e"][0] = dec_min
    er["dec_min_l"][0] = dec_min
    del er["dec_min_l"][1:]
    del er["dec_min_e"][1:]

    er["dec_max_e"][0] = dec_max
    er["dec_max_l"][0] = dec_max
    del er["dec_max_l"][1:]
    del er["dec_max_e"][1:]

    fea, feu, pear = recalculate_metrics(partial, full_400)

    er["full_exposure_area"][0] = fea
    del er["full_exposure_area"][1:]

    er["full_exposure_uncircularity"][0] = feu
    del er["full_exposure_uncircularity"][1:]

    er["partial_exposure_area_ratio"][0] = pear
    del er["partial_exposure_area_ratio"][1:]



def recalculate_bad_bounds(er):
    recalc_eclipse = False

    for i in range(len(er["leg"])):
        if er["ra_min_l"][i] == -180.0 and er["ra_max_l"][i] == 180.0:
            recalc_eclipse = True
            # the other four parameters should be fine as is,
            # and changing as little as possible is good here
            partial = shapely.from_wkb(er["partial"][i])
            if not partial.is_empty:
                _, ra_min, ra_max, _, _, _ = \
                    recalculate_bounds(partial)
            else:
                full = shapely.from_wkb(er["full_400"][i])
                _, ra_min, ra_max, _, _, _ = \
                    recalculate_bounds(full)

            er["ra_min_l"][i] = ra_min
            er["ra_max_l"][i] = ra_max

    recalc_eclipse |= any(emin == -180.0 and emax == 180.0
                          for emin, emax in zip(er["ra_min_e"], er["ra_max_e"]))

    if recalc_eclipse:
        # again, only the ra_min and ra_max parameters are a problem
        partial_e = shapely.union_all([
            shapely.from_wkb(blob) for blob in er["partial"]
        ])
        if not partial_e.is_empty:
            _, ra_min, ra_max, _, _, _ = \
                recalculate_bounds(partial_e)
        else:
            full_e = shapely.intersection_all([
                shapely.from_wkb(blob) for blob in er["full_370"]
            ]).buffer(0)
            _, ra_min, ra_max, _, _, _ = \
                recalculate_bounds(partial_e)

        er["ra_min_e"] = [ra_min]*len(er["ra_min_e"])
        er["ra_max_e"] = [ra_max]*len(er["ra_max_e"])


def make_arrow_table_fragments(er, schema_m, schema_b, schema_l):
    return (
        # only one row per eclipse in the metadata table
        # by construction, all entries in each of these arrays are equal
        Table.from_pydict({
            "eclipse": er["eclipse"][:1],
            "plan_type": er["plan_type"][:1],
            "plan_subtype": er["plan_subtype"][:1],
            "visit": er["visit"][:1],
            "plan_id": er["plan_id"][:1],
            "planned_legs": er["planned_legs"][:1],
            "observed_legs": er["observed_legs"][:1],
            "has_aspect": er["has_aspect"][:1],
            "eclipse_start": er["eclipse_start"][:1],
            "eclipse_duration": er["eclipse_duration"][:1],
            "ok_exposure_time": er["ok_exposure_time"][:1],
            "ra0": er["ra0_e"][:1],
            "ra_min": er["ra_min_e"][:1],
            "ra_max": er["ra_max_e"][:1],
            "dec0": er["dec0_e"][:1],
            "dec_min": er["dec_min_e"][:1],
            "dec_max": er["dec_max_e"][:1],
            "nuv_det_on_time": er["nuv_det_on_time"][:1],
            "nuv_det_temp": er["nuv_det_temp"][:1],
            "nuv_tdc_temp": er["nuv_tdc_temp"][:1],
            "nuv_has_raw6": er["nuv_has_raw6"][:1],
            "fuv_det_on_time": er["fuv_det_on_time"][:1],
            "fuv_det_temp": er["fuv_det_temp"][:1],
            "fuv_tdc_temp": er["fuv_tdc_temp"][:1],
            "fuv_has_raw6": er["fuv_has_raw6"][:1],
        }, schema=schema_m),
        Table.from_pydict({
            "eclipse": er["eclipse"],
            "leg": er["leg"],
            "time": er["time"],
            "duration": er["duration"],
            "ra0": er["ra0_l"],
            "ra_min": er["ra_min_l"],
            "ra_max": er["ra_max_l"],
            "dec0": er["dec0_l"],
            "dec_min": er["dec_min_l"],
            "dec_max": er["dec_max_l"],
            "planned_ra": er["planned_ra"],
            "planned_dec": er["planned_dec"],
            "full_exposure_area": er["full_exposure_area"],
            "full_exposure_uncircularity": er["full_exposure_uncircularity"],
            "partial_exposure_area_ratio": er["partial_exposure_area_ratio"],
        }, schema=schema_b),
        Table.from_pydict({
            "eclipse": er["eclipse"],
            "leg": er["leg"],
            "partial": er["partial"],
            "full_400": er["full_400"],
            "full_370": er["full_370"],
            "full_350": er["full_350"],
        }, schema=schema_l),
    )




def fix_eclipse(
    eclipse_records: dict[str, Any],
    eclipses_to_merge: frozenset[int],
    schema_m, schema_b, schema_l
) -> tuple[Table, Table, Table]:

    eclipse = eclipse_records.pop("eclipse")
    eclipse_records = {
        k.removesuffix("_list"): v
        for k, v in eclipse_records.items()
    }

    if eclipse in eclipses_to_merge:
        merge_legs_and_recalculate_bounds(eclipse_records)
    else:
        recalculate_bad_bounds(eclipse_records)

    return make_arrow_table_fragments(eclipse_records,
                                      schema_m, schema_b, schema_l)

def fix(
    metadata: Table,
    boresight: Table,
    leg_aperture: Table,
    schema_m,
    schema_b,
    schema_l,
    eclipses_to_merge: frozenset[int],
) -> tuple[list[Table], list[Table], list[Table]]:

    mbl = metadata.join(
        boresight,
        keys=["eclipse"],
        left_suffix="_e", right_suffix="_l"
    ).join(
        leg_aperture,
        keys=["eclipse", "leg"],
        left_suffix="_mb", right_suffix="_la"
    )
    check_columns(mbl, EXPECTED_COLS_AFTER_JOINS)

    new_metadata = []
    new_boresight = []
    new_leg_aperture = []

    # pyarrow should have "for group in table.groupby()"...
    # pandas does but I want to avoid converting to pandas boxes
    # (we can live with conversion to native python objects)
    for eclipse in mbl.group_by(["eclipse"]).aggregate(
        [(name, "list") for name in mbl.column_names]
    ).to_pylist():
        nm, nb, nl = fix_eclipse(eclipse, eclipses_to_merge,
                                 schema_m, schema_b, schema_l)
        new_metadata.append(nm)
        new_boresight.append(nb)
        new_leg_aperture.append(nl)


    return (
        new_metadata,
        new_boresight,
        new_leg_aperture,
    )


def vstack_and_sort(
    tables: list[Table],
    sort_order: list[tuple[str, str]]
) -> Table:
    return concat_tables(tables).sort_by(sort_order)


def rewrite_table(path: Path, table: Table, *args: Any, **kwargs: Any) -> None:
    """Replace the file named "path" with a parquet file containing
       the contents of "table".  If "path" doesn't exist, create it,
       otherwise rename it with a ".bak" suffix; in that case, any
       existing file named "{path}.bak" will be overwritten."""
    target_dir = path.parent
    target_stem = path.stem

    backup_path = path.with_name(path.name + ".bak")

    # Infuriatingly, there is no way to create an Arrow NativeFile
    # directly from a file descriptor.  (You can create one from a
    # Python file object, but then you have multiple layers of buffering.)
    tempfd, temppath = tempfile.mkstemp(
        dir=target_dir,
        prefix=target_stem,
        suffix=".partial",
    )
    os.close(tempfd)
    temppath = Path(temppath)
    try:
        parquet.write_table(table, temppath, *args, **kwargs)
        try:
            path.rename(backup_path)
        except FileNotFoundError:
            pass

        temppath.rename(path)

    except:
        try:
            temppath.unlink()
        except FileNotFoundError:
            pass
        raise


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("datadir", type=Path)
    ap.add_argument("eclipse_to_merge", type=int, nargs="*")
    ap.add_argument("--only-fixed", action="store_true")

    args = ap.parse_args()
    eclipses_to_merge = frozenset(args.eclipse_to_merge)
    metadata_path     = args.datadir / "metadata.parquet"
    boresight_path    = args.datadir / "boresight.parquet"
    leg_aperture_path = args.datadir / "leg-aperture.parquet"

    filter_boresight_to_fix = (
        pc.field("eclipse").isin(eclipses_to_merge)
        | ((pc.field("ra_min") == -180.0) & (pc.field("ra_max") == 180.0))
    )

    boresight_to_fix = parquet.read_table(
        boresight_path, filters = filter_boresight_to_fix
    )
    boresight_okay = parquet.read_table(
        boresight_path, filters = ~filter_boresight_to_fix
    )
    check_columns(boresight_to_fix, EXPECTED_COLS_BORESIGHT)
    check_columns(boresight_okay, EXPECTED_COLS_BORESIGHT)

    schema_b = boresight_to_fix.schema
    assert boresight_okay.schema == schema_b

    filter_eclipses_to_fix = pc.field("eclipse").isin(frozenset(
        boresight_to_fix.column("eclipse").to_pylist()
    ))

    metadata_to_fix = parquet.read_table(
        metadata_path, filters = filter_eclipses_to_fix
    )
    metadata_okay = parquet.read_table(
        metadata_path, filters = ~filter_eclipses_to_fix
    )
    check_columns(metadata_to_fix, EXPECTED_COLS_METADATA)
    check_columns(metadata_okay, EXPECTED_COLS_METADATA)
    schema_m = metadata_to_fix.schema
    assert metadata_okay.schema == schema_m

    leg_aperture_to_fix = parquet.read_table(
        leg_aperture_path, filters = filter_eclipses_to_fix
    )
    leg_aperture_okay = parquet.read_table(
        leg_aperture_path, filters = ~filter_eclipses_to_fix
    )
    check_columns(leg_aperture_to_fix, EXPECTED_COLS_LEG_APERTURE)
    check_columns(leg_aperture_okay, EXPECTED_COLS_LEG_APERTURE)
    schema_l = leg_aperture_to_fix.schema
    assert leg_aperture_okay.schema == schema_l

    metadata_fixed, boresight_fixed, leg_aperture_fixed = \
        fix(metadata_to_fix, boresight_to_fix, leg_aperture_to_fix,
            schema_m, schema_b, schema_l,
            eclipses_to_merge)

    if args.only_fixed:
        metadata_path     = args.datadir / "metadata-fixed.parquet"
        boresight_path    = args.datadir / "boresight-fixed.parquet"
        leg_aperture_path = args.datadir / "leg-aperture-fixed.parquet"
    else:
        metadata_fixed.append(metadata_okay)
        boresight_fixed.append(boresight_okay)
        leg_aperture_fixed.append(leg_aperture_okay)

    sort_order_m = [("eclipse", "ascending")]
    sort_order_b = [("eclipse", "ascending"), ("leg", "ascending")]
    metadata_new = vstack_and_sort(metadata_fixed, sort_order_m)
    boresight_new = vstack_and_sort(boresight_fixed, sort_order_b)
    leg_aperture_new = vstack_and_sort(leg_aperture_fixed, sort_order_b)

    rewrite_table(metadata_path, metadata_new,
                  write_page_checksum=True,
                  sorting_columns=SC.from_ordering(
                      metadata_new.schema,
                      sort_order_m
                  ))

    rewrite_table(boresight_path, boresight_new,
                  write_page_checksum=True,
                  sorting_columns=SC.from_ordering(
                      boresight_new.schema,
                      sort_order_b
                  ))

    # this one has really big blobs and is known not to benefit from
    # dictionary encoding
    rewrite_table(leg_aperture_path, leg_aperture_new,
                  write_page_checksum=True,
                  use_dictionary=False, compression="zstd",
                  sorting_columns=SC.from_ordering(
                      leg_aperture_new.schema,
                      sort_order_b
                  ))


main()
