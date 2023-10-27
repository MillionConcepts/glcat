import sys

from dustgoggles.dynamic import exc_report
from dustgoggles.structures import MaybePool

sys.path.insert(0, '/home/bekah/gPhoton2')
from gPhoton.types import Pathlike, GalexBand
from typing import Optional, Literal
from pyarrow import parquet
import pandas as pd
from aspect_correction.slew_correction import refine_slew_frame
from aspect_correction.dose_aspect_correction import refine_normal_frame
from gPhoton.reference import PipeContext
from backplanes import dosemaps_just_for_timestamps
# aspect root = '/home/bekah/gPhoton2'/test_data'
# gphoton root = '/home/bekah/gPhoton2'


def refine_eclipse(
        eclipse: int,
        xylist: bool,
        band: Literal["NUV", "FUV"],
        aspect_root: str,
        gphoton_root: str,
        ext="gzip",
        threads=None,
        tmp_root="/tmp"):
    """ main pipeline for processing an eclipse frame by frame.
     this pipeline requires that all 1s backplanes be pre-produced.
      easiest script would run backplanes for a particular eclipse at 1s
       intervals and then runs this function.
       :param eclipse: eclipse number, not paddded
       :param xylist: true if wantt to use xylists for norm frame processing
       :param band: 'NUV' or 'FUV'
       :param aspect_root: where backplanes are and you want astrometry to save
       output
       :param gphoton_root: where gphoton is, needed for metadata files etc
       :param ext: file compression
       :param threads: number of frames to run in parallel; None for serial
       :param tmp_root: where temp files are put by astrom, need for cleanup"""

    # setup to run pipeline by getting relevant info & paths
    metadata_paths = metadata_filepaths(gphoton_root)
    eclipse_info = get_eclipse_info(eclipse, metadata_paths, band)

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
    eclipse_aspect = {}

    if eclipse_info['actual_legs'] == 0:
        modified_frame_list = get_backplane_filenames(
            eclipse_info,
            frame_list,
            paths[0],
            0,
            aspect_root)
    for leg in range(eclipse_info['actual_legs']):
        files = get_backplane_filenames(
            eclipse_info,
            frame_list,
            paths[leg],
            leg,
            aspect_root)
        modified_frame_list = pd.concat([modified_frame_list, files], axis=0)
    pool = MaybePool(threads)
    base_args = [eclipse_info, aspect_root, xylist]
    frame_args = [
        {'args': [modified_frame_list.iloc[frame]] + base_args}
        for frame in range(len(modified_frame_list))
    ]
    # print(f"Running refine on frame {frame}")
    # NOTE: will want to put statements like this inside refine_frame
    try:
        pool.map(refine_frame, frame_args)
        pool.close()
        pool.join()
        aspects = pool.get()
    finally:
        pool.terminate()
    for frame, aspect in aspects.items():
        try:
            if isinstance(aspect, Exception):
                raise aspect
            elif aspect is None:
                raise ValueError("aspect is None")
            eclipse_aspect[frame] = {
                'ra': float(aspect.iloc[0]['ra_tangent']),
                'dec': float(aspect.iloc[0]['dec_tangent']),
                'roll': float(aspect.iloc[0]['orientation']),
                'pixscale': float(aspect.iloc[0]['pixscale']),
                'ra_center': float(aspect.iloc[0]['ra_center']),
                'dec_center': float(aspect.iloc[0]['dec_center']),
                'logodds': float(aspect.iloc[0]['logodds'])
            }
        except KeyboardInterrupt:
            raise
        except Exception as ex:
            print(f"Error in frame {frame}")
    # convert dict of aspect solns to pd df and return
    aspect_soln = pd.DataFrame.from_dict(eclipse_aspect, orient='index')
    aspect_soln.to_csv(f"{aspect_root}/e{eclipse_info['eclipse_str']}"
                       f"/{eclipse_info['eclipse_str']}_new_aspect.csv")
    # clean up
    cleanup_eclipse_files(aspect_root, eclipse_info)
    cleanup_temp_files(tmp_root)

    return aspect_soln


def cleanup_eclipse_files(aspect_root, eclipse_info):
    """ delete aspect directory (not final asp soln) and also the intermediate
    files created by astrometry.net """
    import shutil
    folder_path = f"{aspect_root}+{eclipse_info['eclipse_str']}/aspect"
    try:
        shutil.rmtree(folder_path)
        print(f"Deleted folder: {folder_path}")
    except OSError as e:
        print(f"Error deleting folder: {e}")
    return


def cleanup_temp_files(tmp_root):
    """ remove tmp files produced by astrometry.net """
    import os
    strings_to_check = ["tmp.uncompressed", "tmp.ppm", ".pnm", "tmp.axy"]
    for filename in os.listdir(tmp_root):
        file_path = os.path.join(tmp_root, filename)
        if any(substring in filename for substring in strings_to_check):
            try:
                os.remove(file_path)
                print(f"Removed file: {file_path}")
            except OSError as e:
                print(f"Error removing file: {e}")
    return


def refine_frame(frame_series, eclipse_info, aspect_root, xylist):
    """refines a frame based on info in frame series (ex: normal or slew frame).
    important cols include: time, flags_x,ptag,hvnom_nuv,hvnom_fuv,ra_acs,
    dec_acs, roll_acs,frame_type,backplane_path,time_stamp,leg """
    if frame_series['frame_type'] == "slew":
        print(f"Running slew frame {frame_series['time']}.")
        aspect = refine_slew_frame(
            frame_series,
            eclipse_info,
            aspect_root)
    elif frame_series['frame_type'] == "ref":
        print(f"Running ref frame {frame_series['time']}.")
        aspect = refine_normal_frame(frame_series, xylist)
    else:
        print(f"Frame type not accepted, frame {frame_series['time']}.")
        aspect = ()
    aspect['time'] = frame_series['time']
    aspect['frame_type'] = frame_series['frame_type']
    aspect['flags_x'] = frame_series['flags_x']
    aspect['hvnom_nuv'] = frame_series['hvnom_nuv']
    aspect['hvnom_fuv'] = frame_series['hvnom_fuv']
    aspect['mission_time'] = frame_series['mission_time']
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
        metadata_paths: dict,
        band: Literal["NUV", "FUV"]):
    """returns a dictionary of info about eclipse"""
    from gPhoton.reference import titular_legs
    actual, nominal = titular_legs(eclipse)
    eclipse_info = {'actual_legs': actual,
                    'nominal_legs': nominal,
                    'band': band}
    metadata = parquet.read_table(metadata_paths['metadata'],
                                  filters=[('eclipse', '==', eclipse)]).to_pandas()
    eclipse_info['ra_max'] = metadata['ra_max']
    eclipse_info['ra_min'] = metadata['ra_min']
    eclipse_info['dec_max'] = metadata['dec_max']
    eclipse_info['dec_min'] = metadata['dec_min']
    eclipse_info['obstype'] = metadata['obstype']
    eclipse_info['eclipse'] = eclipse
    eclipse_info['eclipse_str'] = str(eclipse).zfill(5)
    return eclipse_info


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
    from gPhoton.reference import eclipse_to_paths
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


def get_frame_list_mod(
        eclipse: int,
        metadata_paths: dict):
    """ join extended and og aspect parquet per eclipse to get unrefined
    'slew' frames. MODIFIED VERSION to use parts of backplane code """
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
        leg: int,
        aspect_root: str):
    """ produce df of backplane filenames + other relevant names that are
     frame specific (backplane, xylist, wcs etc) """
    phot = parquet.read_table(
        paths['photonfile'], columns=['t', 'col', 'row', 'detrad']
    )
    # select frames that make up a leg
    t_f = phot['t'][0].as_py()  # first timestamp of leg from photonlist
    t_l = phot['t'][-1].as_py()  # last timestamp of leg from photonlist

    # rounding to have file names match backplane names (not sure why
    # but timestamps are less sig figs in the aspect table)
    t_f = truncate(t_f)
    t_l = truncate(t_l)
    #frame_list = frame_list.round({'time': 3})
    frame_list['time'] = frame_list['time'].apply(truncate)

    leg_frames = frame_list.loc[(frame_list['time'] >= t_f) & (frame_list['time'] <= t_l)].copy()

    # backplanes name formatting fNNNNdd_tNNNNdd
    # this is probably an overly complex way to do the naming but ... it works
    leg_frames['backplane_path_og'] = paths['movie'].replace('.fits.gz', f'_dose.fits') \
        .replace('movie', 'tmovie')
    # subtract first timestamp of photonlist
    leg_frames['mission_time'] = leg_frames['time'].copy()
    leg_frames['time'] = leg_frames['time'] - t_f
    leg_frames['time_stamp'] = leg_frames['time'].astype(str).str.split('.').str[0] \
        .str.zfill(4).str.cat(leg_frames['time'].astype(str).str.split('.').str[1].str[:4])
    # backplane path
    leg_frames['backplane_path'] = leg_frames['backplane_path_og'].str.split('movie').str[0] \
        .str.cat(leg_frames['time_stamp']).str.cat(leg_frames['backplane_path_og']
                                                   .str.split('movie').str[1]) + ".gz"
    # leg
    leg_frames['leg'] = leg
    # xylist path = l[leg]ts[time stamp].xyls
    leg_frames['aspect_output'] = aspect_root + f"/e{eclipse_info['eclipse_str']}/astrom"
    leg_frames['xylist_path'] = leg_frames['aspect_output'] + "/l" + leg_frames['leg'].astype(str)\
                            + "ts" + (leg_frames['time_stamp'].astype(str)) + ".xyls"
    # wcs path
    leg_frames['wcs_path'] = leg_frames['aspect_output'] + "/l" + leg_frames['leg'].astype(str)\
                            + "ts" + (leg_frames['time_stamp'].astype(str)) + ".wcs"
    return leg_frames


def get_backplane_filenames_mod(
        eclipse_info: dict,
        frame_list,
        paths: dict,
        leg: int,
        aspect_root: str):
    """ produce df of backplane filenames + other relevant names that are
     frame specific (backplane, xylist, wcs etc) """
    eclipse = eclipse_info["eclipse"]
    band = eclipse_info["band"]
    depth = 1
    leg = leg
    threads = 4
    burst = True
    local = "/media/bekah/BekahA/glcat"
    kind = "dose"
    radius = 600
    write = {'array': True, 'xylist': False}
    inline = True
    threshold = .45
    star_size = 2
    write = {} if write is None else dict(write)
    ctx = PipeContext(
        eclipse,
        band,
        depth,
        "gzip",
        local,
        leg=leg,
        threads=threads,
        burst=burst,
        write=write,
        start_time=1000,
        stop_after=None,
    )
    ctx.hdu_constructor_kwargs = dict(ctx.hdu_constructor_kwargs)
    ctx.inline = inline
    ctx.star_size = star_size
    ctx.threshold = threshold

    # merge backplane list with naming df
    time_stamps = dosemaps_just_for_timestamps(ctx, radius=400)
    leg_frames = frame_list.merge(time_stamps, on="time", how="outer")

    # backplanes name formatting fNNNNdd_tNNNNdd
    # this is probably an overly complex way to do the naming but ... it works
    leg_frames['backplane_path_og'] = paths['movie'].replace('.fits.gz', f'_dose.fits') \
        .replace('movie', 'tmovie')
    # subtract first timestamp of photonlist
    leg_frames['time_stamp'] = leg_frames['mod_time'].astype(str).str.split('.').str[0] \
        .str.zfill(4).str.cat(leg_frames['mod_time'].astype(str).str.split('.').str[1].str[:4])
    # backplane path
    leg_frames['backplane_path'] = leg_frames['backplane_path_og'].str.split('movie').str[0] \
        .str.cat(leg_frames['time_stamp']).str.cat(leg_frames['backplane_path_og']
                                                   .str.split('movie').str[1])
    # leg
    leg_frames['leg'] = leg
    # xylist path = l[leg]ts[time stamp].xyls
    leg_frames['aspect_output'] = aspect_root + f"/e{eclipse_info['eclipse_str']}/astrom"
    leg_frames['xylist_path'] = leg_frames['aspect_output'] + "/l" + leg_frames['leg'].astype(str)\
                            + "ts" + (leg_frames['time_stamp'].astype(str)) + ".xyls"
    # wcs path
    leg_frames['wcs_path'] = leg_frames['aspect_output'] + "/l" + leg_frames['leg'].astype(str)\
                            + "ts" + (leg_frames['time_stamp'].astype(str)) + ".wcs"
    return leg_frames


def truncate(num):
    """ literally just for cutting floats to 3
    dec places without rounding. """
    places = 3
    mult = 10 ** places
    return int(num * mult) / mult

