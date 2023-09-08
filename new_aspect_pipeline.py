import sys
sys.path.insert(0, '/home/bekah/gPhoton2')
from gPhoton.reference import get_legs, titular_legs, check_eclipse, eclipse_to_paths
from gPhoton.types import Pathlike, GalexBand
from typing import Any, Callable, Sequence, Mapping, Optional, Literal
from pyarrow import parquet
import pandas as pd
from aspect_correction.dose_aspect_correction import refine_normal_frame
from slew_aspect_correction.main import refine_slew_frame


# aspect root = '/home/bekah/gPhoton2'/test_data'
# gphoton root = '/home/bekah/gPhoton2'


def refine_eclipse(
        eclipse: int,
        band: str,
        aspect_root: str,
        gphoton_root: str,
        ext="gzip"):
    """ main pipeline for processing an eclipse frame by frame """

    # setup to run pipeline by getting relevant info & paths
    metadata_paths = metadata_filepaths(gphoton_root)
    eclipse_info = get_eclipse_info(eclipse, metadata_paths)

    # paths for photonlist, possible output images etc
    # function is from gPhoton2
    paths = get_base_file_paths(
        eclipse=eclipse,
        band=band,
        ext=ext,
        root=aspect_root,
        mode='direct',
        legs=eclipse_info['actual_legs'])

    # get list of frames for all legs into a single dataframe,
    # with paths for each backplane
    frame_list = get_frame_list(eclipse, metadata_paths)
    modified_frame_list = pd.DataFrame()
    if eclipse_info['actual_legs'] == 0:
        modified_frame_list = get_backplane_filenames(
            eclipse_info,
            frame_list,
            paths[0],
            0)
    for leg in range(eclipse_info['actual_legs']):
        files = get_backplane_filenames(
            eclipse_info,
            frame_list,
            paths[leg],
            leg)
        modified_frame_list = pd.concat([modified_frame_list, files], axis=0)

    for frame in range(len(modified_frame_list)):
        aspect = refine_frame(modified_frame_list.iloc[frame])
        # aspect is tuple of ra, dec, roll, time
        3  # TODO: do I want the new aspect to go in the old mod_frame_list or it's own file?
        # modified_frame_list.iloc
    return


def refine_frame(frame_series):
    """refines a frame based on info in frame series (ex: normal or slew frame).
    important cols include: time, flags_x,ptag,hvnom_nuv,hvnom_fuv,ra_acs,
    dec_acs, roll_acs,frame_type,backplane_path,time_stamp,leg """
    if frame_series['frame_type'] == "slew":
        aspect = refine_slew_frame(frame_series)
    elif frame_series['frame_type'] == "ref":
        aspect = refine_normal_frame(frame_series)
    return aspect


def metadata_filepaths(root: str):
    """return paths to aspect, boresight, and metadata parquet files """
    metadata_paths = {'og_aspect': root + '/gPhoton/aspect/aspect.parquet',
                      'expanded_aspect': root + '/gPhoton/aspect/aspect2.parquet',
                      'astrometry_aspect': root + '/gPhoton/aspect/astrom_aspect.parquet',
                      'boresight': root + '/gPhoton/aspect/boresight.parquet',
                      'metadata': root + '/gPhoton/aspect/metadata.parquet'}
    # astrom aspect doesn't exist yet but might use it
    return metadata_paths


def get_eclipse_info(
        eclipse: int,
        metadata_paths: dict):
    """returns a dictionary of info about eclipse"""
    actual, nominal = titular_legs(eclipse)
    eclipse_info = {'actual_legs': actual, 'nominal_legs': nominal}
    metadata = parquet.read_table(metadata_paths['metadata'],
                                  filters=[('eclipse', '==', eclipse)]).to_pandas()
    eclipse_info['ra_max'] = metadata['ra_max']
    eclipse_info['ra_min'] = metadata['ra_min']
    eclipse_info['dec_max'] = metadata['dec_max']
    eclipse_info['dec_min'] = metadata['dec_min']
    eclipse_info['obstype'] = metadata['obstype']
    return eclipse_info


# TODO: need to check if legs naming starts with 0 or 1
def get_base_file_paths(
        eclipse: int,
        band: GalexBand = "NUV",
        depth: Optional[int] = None,
        compression: Literal["none", "gzip", "rice"] = "gzip",
        root: Pathlike = "data",
        start: Optional[float] = None,
        mode: str = "direct",
        legs: int = 0,
        aperture: Optional[float] = None,
        **kwargs,
) -> dict[str, str]:
    """dictionary of dictionaries of file paths for each leg of eclipse """
    paths = {}
    if legs == 0:
        paths[0] = eclipse_to_paths(eclipse=eclipse, band=band, depth=1,
                                    compression=compression, root=root,
                                    mode=mode, leg=legs, )
    else:
        for i in range(legs):
            paths[i] = eclipse_to_paths(eclipse=eclipse, band=band, depth=1,
                                        compression=compression, root=root,
                                        mode=mode, leg=i)
    return paths


def get_frame_list(
        eclipse: int,
        metadata_paths: dict):
    """ join extended and og aspect parquet per eclipse to get unrefined
    'slew' frames """
    from pyarrow import parquet
    og_aspect = parquet.read_table(metadata_paths['og_aspect'],
                                   filters=[('eclipse', '==', eclipse)]).to_pandas()
    og_aspect['original'] = 'ref'

    new_aspect = parquet.read_table(metadata_paths['expanded_aspect'],
                                    filters=[('eclipse', '==', eclipse)]).to_pandas()
    new_aspect = new_aspect.rename(columns={"pktime": "time"})
    new_aspect['original'] = 'slew'

    new_aspect = new_aspect.merge(og_aspect, on="time", how="outer")

    new_aspect['original_y'].fillna('slew', inplace=True)
    new_aspect = new_aspect.rename(columns={"original_y": "frame_type"})

    # get rid of weird distant time stamp at end of file sometimes
    # (only works if it's one entry)
    if new_aspect.iloc[-1].time - new_aspect.iloc[-2].time > 10:
        new_aspect.drop(new_aspect.tail(1).index, inplace=True)
    return new_aspect


def get_backplane_filenames(
        eclipse_info: dict,
        frame_list,
        paths: dict,
        leg: int):
    """ produce df of backplane filenames """
    phot = parquet.read_table(
        paths['photonfile'], columns=['t', 'col', 'row', 'detrad']
    )
    # select frames that make up a leg
    t_f = phot['t'][0].as_py()  # first timestamp of leg from photonlist
    t_l = phot['t'][-1].as_py()  # last timestamp of leg from photonlist
    leg_frames = frame_list.loc[(frame_list['time'] >= t_f) & (frame_list['time'] <= t_l)].copy()
    # backplanes name formatting fNNNNdd_tNNNNdd
    # this is probably an overly complex way to do the naming but ... it works
    leg_frames['backplane_path'] = paths['movie'].replace('.fits.gz', f'_dose.fits') \
        .replace('movie', 'tmovie')
    # subtract first timestamp of photonlist
    leg_frames['time'] = leg_frames['time'] - t_f
    leg_frames['time_stamp'] = leg_frames['time'].astype(str).str.split('.').str[0] \
        .str.zfill(4).str.cat(leg_frames['time'].astype(str).str.split('.').str[1].str[:4])
    # backplane path
    leg_frames['backplane_path'] = leg_frames['backplane_path'].str.split('movie').str[0] \
        .str.cat(leg_frames['time_stamp']).str.cat(leg_frames['backplane_path']
                                                   .str.split('movie').str[1])
    # leg
    leg_frames['leg'] = leg
    # xylist path = l[leg]ts[time stamp].xyls
    leg_frames['xylist_path'] = "l" + leg_frames['leg'].astype(str) + "ts" + \
                                (leg_frames['time_stamp'].astype(str)) + ".xyls"

    return leg_frames


