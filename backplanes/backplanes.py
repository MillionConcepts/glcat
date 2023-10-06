from functools import partial
from itertools import product
from multiprocessing import Pool
from pathlib import Path
from types import MappingProxyType
from typing import Literal, Optional, Mapping

import astropy.io.fits
import fast_histogram as fh
import numpy as np
import pandas as pd
import sh
import sys
#sys.path.append('/home/bekah/gphoton_working')
sys.path.append('/home/bekah/gPhoton2')
sys.path.append('/home/ubuntu/gPhoton2')

from gPhoton.types import GalexBand
from more_itertools import windowed
from pyarrow import parquet
from scipy.stats import binned_statistic_2d
from gPhoton.moviemaker._steps import (
    select_on_detector,
    generate_wcs_components,
    slice_into_memory, populate_fits_header
)
from gPhoton.pretty import print_inline
from gPhoton.reference import PipeContext
from gPhoton.sharing import (reference_shared_memory_arrays)
from gPhoton.vorpal import between
from aspect_correction.dose_aspect_correction import get_stars


def stub_header(band, wcs=None, tranges=None):
    """
    create an astropy.io.fits.Header object containing our canonical
    metadata values
    """
    header = astropy.io.fits.Header()
    if wcs is not None:
        header["CDELT1"], header["CDELT2"] = wcs.wcs.cdelt
        header["CTYPE1"], header["CTYPE2"] = wcs.wcs.ctype
        header["CRPIX1"], header["CRPIX2"] = wcs.wcs.crpix
        header["CRVAL1"], header["CRVAL2"] = wcs.wcs.crval
        header["EQUINOX"], header["EPOCH"] = 2000.0, 2000.0
    header["BAND"] = 1 if band == "NUV" else 2
    if tranges is not None:
        for i, trange in enumerate(tranges):
            header["EXPSTART"] = np.array(tranges).min()
            header["EXPEND"] = np.array(tranges).max()
            header["N_FRAME"] = len(tranges)
            header["T0_{i}".format(i=i)] = trange[0]
            header["T1_{i}".format(i=i)] = trange[1]
    return header


def send_to_shared_memory(components, depth):
    total_trange = (components['t'].min(), components['t'].max())
    t0s = np.arange(total_trange[0], total_trange[1] + depth, depth)
    tranges = list(windowed(t0s, 2))
    ax_blocks = {}
    for frame_ix, trange in enumerate(tranges):
        frame_time_ix = between(components['t'], *tranges[frame_ix])[0]
        ax_blocks[frame_ix] = slice_into_memory(
            {k: v for k, v in components.items() if k != 't'},
            (frame_time_ix.min(), frame_time_ix.max())
        )
    return ax_blocks, tranges


def load_for_dosemap(photonfile, radius=400, snippet: Optional[tuple] = None):
    """ modified to allow for specific t-range selection (snippet) within
     parquet file (for slew frames, for example) """
    if snippet is None:
        phot = parquet.read_table(
            photonfile, columns=['t', 'col', 'row', 'detrad']
        )
        t = phot['t'][0]
    else:
        # filter to load just rows for a specific time range
        phot = parquet.read_table(
            photonfile, columns=['t', 'col', 'row', 'detrad']
        )
        t = phot['t'][0]
        phot = parquet.read_table(
            photonfile, columns=['t', 'col', 'row', 'detrad'],
            filters=[('t', '>', snippet[0]-1), ('t', '<', snippet[1]+1)]
        )
    phot = select_on_detector(phot, radius)
    return {
        't': phot['t'].to_numpy(),
        'x': phot['row'].to_numpy() * 4,
        'y': phot['col'].to_numpy() * 4
    }, t


def dosemap_ranges(radius):
    slop = 400 - radius
    return [slop, 3200 - slop], [slop, 3200 - slop]


def dosemap_frame(axes, radius=400, resolution=None):
    """
    make a 'dosemap' from a photonfile based on col/row.
    the photonfile must have been written in 'extended' mode.
    no photons outside radius will be rendered.
    400 should in general capture all 'real' photons.
    750 should usually be safe to capture the stims as well.
    if resolution is None, set the resolution (X x Y)
    of the output array equal to 8 * radius.
    axes of returned array are col/row upscaled by 4.
    """
    resolution = radius * 8 if resolution is None else radius
    return fh.histogram2d(
        axes['x'], axes['y'], bins=resolution, range=dosemap_ranges(radius)
    ).astype('float32')


def sm_make_dosemap(block_info, radius, frame_ix):
    _, axes = reference_shared_memory_arrays(block_info)
    print(f'integrating frame {frame_ix} dosemap')
    return {'dose': dosemap_frame(axes, radius)}


def sm_compute_dosemap_frame(block_info, imsz, frame_ix, ctx, start_time, time_stamp):
    maps = sm_make_dosemap(block_info, imsz, frame_ix)
    for _, block in reference_shared_memory_arrays(
        block_info, fetch=False
    )[0].items():
        block.close()
        block.unlink()
    if ctx.write.get('xylist') is True:
        write_xylist_inline(ctx, frame_ix, maps)
    return write_or_return_arrays(maps, frame_ix, ctx, start_time, time_stamp)


def write_xylist_inline(ctx, frame_ix, maps):
    # TODO: let's unravel this later.
    file_names = {
        'xylist': {
            frame_ix: str(Path(ctx.eclipse_path(), f'frame{frame_ix}.xyls'))
        }
    }
    table_name = get_stars(
        maps['dose'], frame_ix, file_names, ctx.threshold, ctx.star_size
    )
    print(f"wrote xylist for frame {frame_ix} to {table_name}")


def make_dosemap(ctx: PipeContext, radius: int = 400):
    components, start_time = load_for_dosemap(ctx()['photonfile'], radius, ctx.snippet)
    if ctx.depth is None:
        maps = dosemap_frame(components, radius)
        if ctx.write['array'] is True:
            if ctx.write['xylist'] is True:
                print('note: xylists not implemented for full-depth dosemaps.')
            return write_backplane_image(maps, ctx, start_time)
    print('slicing position data into shared memory')
    dose_blocks, tranges = send_to_shared_memory(components, ctx.depth)
    del components
    frames = {}
    time_stamps = {}
    pool = Pool(ctx.threads) if ctx.threads is not None else None
    for frame_ix, trange in enumerate(tranges):
        # frame value is trange[0] - start time of eclipse
        time_stamps[frame_ix] = trange[0]-start_time.as_py()
        if pool is not None:
            frames[frame_ix] = pool.apply_async(
                sm_compute_dosemap_frame,
                (dose_blocks[frame_ix], radius, frame_ix, ctx, start_time, time_stamps[frame_ix])
            )
        else:
            frames[frame_ix] = sm_compute_dosemap_frame(
                dose_blocks[frame_ix], radius, frame_ix, ctx, start_time, time_stamps[frame_ix]
            )
    if pool is not None:
        pool.close()
        pool.join()
        frames = {ix: frame.get() for ix, frame in frames.items()}
    ranges = dosemap_ranges(radius)
    imsz = [ranges[0][1] - ranges[0][0]] * 2
    return write_backplane_movies(frames, imsz, ctx, tranges, start_time, time_stamps)


def load_for_xymap(photonfile, radius=400):
    xy_cols = ["t", "ra", "dec", "col", "row", "detrad"]
    xytab = select_on_detector(
        parquet.read_table(photonfile, columns=xy_cols), radius
    )
    foc, wcs = generate_wcs_components(xytab)
    imsz = (
        int((wcs.wcs.crpix[1] - 0.5) * 2), int((wcs.wcs.crpix[0] - 0.5) * 2)
    )
    components = {
        't': xytab['t'].to_numpy(),
        'x': xytab['row'].to_numpy() * 4,
        'y': xytab['col'].to_numpy() * 4,
        'foc': foc
    }
    return components, wcs, imsz


def make_full_depth_xymaps(components, imsz):
    kwargs = {
        'x': components['foc'][:, 1] - 0.5,
        'y': components['foc'][:, 0] - 0.5,
        'bins': imsz,
        'range': ([[0, imsz[0]], [0, imsz[1]]])
    }
    binner = partial(binned_statistic_2d, **kwargs)
    xymaps = {}
    for ax, stat in product(('x', 'y'), ('mean', 'std')):
        print(f'making {ax}_{stat}')
        xymaps[f'{ax}_{stat}'] = binner(
            values=components[ax], statistic=stat
        )[0]
    return xymaps


# noinspection PyTypeChecker
def sm_binner(arrays, ax, stat, imsz):
    return binned_statistic_2d(
        arrays['foc'][:, 1] - 0.5,
        arrays['foc'][:, 0] - 0.5,
        arrays[ax],
        bins=imsz,
        range=([[0, imsz[0]], [0, imsz[1]]]),
        statistic=stat
    )[0].astype(np.float32)


def sm_make_xymaps(block_info, imsz, frame_ix):
    _, arrays = reference_shared_memory_arrays(block_info)
    xymaps = {}
    for ax, stat in product(('x', 'y'), ('mean', 'std')):
        print_inline(f'integrating {frame_ix} {ax}_{stat}')
        xymaps[f'{ax}_{stat}'] = sm_binner(arrays, ax, stat, imsz)
    return xymaps


def sm_compute_xymap_frame(block_info, imsz, frame_ix, ctx, wcs, tranges):
    maps = sm_make_xymaps(block_info, imsz, frame_ix)
    for _, block in reference_shared_memory_arrays(
        block_info, fetch=False
    )[0].items():
        block.close()
        block.unlink()
    return write_or_return_arrays(maps, frame_ix, ctx, wcs, tranges)


def write_backplane_file(
    image, ctx, start_time, tranges=None, wcs=None, name="", frame="movie", time_stamp=0):
    if frame == 'image':
        fn = ctx()['image']
    else:
        #fNNNNdd_tNNNNdd
        fn = ctx(frame=frame)['movie']
        split_time = str(time_stamp).split('.')
        N = split_time[0].zfill(4)
        # first 4 of decimal place
        d = split_time[1][:4]
        fn = fn.replace('movie', f't{N}{d}')
    fn = fn.replace('.fits.gz', f'_{name}.fits')
    for ext in ('', '.gz'):
        if Path(ctx.eclipse_path(), fn + ext).exists():
            Path(fn + ext).unlink()
    header = stub_header(ctx.band, wcs, tranges)
    hdu = astropy.io.fits.PrimaryHDU(image, header=header)
    hdu.writeto(fn)
    print(f'gzipping {name} {frame}')
    sh.igzip(fn)
    Path(fn).unlink()


def write_backplane_image(maps, ctx, start_time, tranges=None, wcs=None):
    for name, image in maps.items():
        write_backplane_file(image, ctx, start_time, tranges, wcs, name, 'image')


def sparse_to_movie(sparse, imsz):
    return np.dstack([f.to_dense().reshape(imsz) for f in sparse])


def write_backplane_movies(frames, imsz, ctx, tranges, wcs=None, start_time=0, time_stamps=None):
    if check_inline_write(ctx) is True and ctx.snippet is not None:
        print("all outputs written inline; returning time ranges for snippet.")
        return dosemaps_just_for_timestamps(ctx, radius=400)
    if check_inline_write(ctx) is True:
        print("all outputs written inline; terminating.")
        return
    if ctx.write.get('array') is False:
        print('write["array"] = False; terminating.')
        return
    args = (ctx, start_time, tranges, wcs)
    if ctx.burst is False:
        maps = {k: [] for k in frames[0].keys()}
        for ix in list(frames.keys()):
            for name, array in frames[ix].items():
                maps[name].append(array)
            del frames[ix]
        for name in list(maps.keys()):
            print(f'writing {name}')
            write_backplane_file(
                sparse_to_movie(maps[name], imsz), *args, name, "movie", time_stamps[name]
            )
            del maps[name]
        return
    for ix in list(frames.keys()):
        for name, sparse in frames[ix].items():
            print(f'writing {name} (frame {ix})')
            write_backplane_file(
                sparse.to_dense().reshape(imsz), *args, name, str(ix).zfill(5), time_stamps[name]
            )
        del frames[ix]


def make_xymaps(ctx: PipeContext):
    print("loading photonlist and computing wcs")
    components, wcs, imsz, start_time = load_for_xymap(ctx()['photonfile'])
    if ctx.depth is None:
        maps = make_full_depth_xymaps(components, imsz)
        return write_backplane_image(maps, ctx, start_time, None, wcs)
    print('slicing position data into shared memory')
    ax_blocks, tranges = send_to_shared_memory(components, ctx.depth)
    del components
    frames = {}
    pool = Pool(ctx.threads) if ctx.threads is not None else None
    for frame_ix, trange in enumerate(tranges):
        if pool is not None:
            frames[frame_ix] = pool.apply_async(
                sm_compute_xymap_frame,
                (ax_blocks[frame_ix], imsz, frame_ix, ctx, wcs, tranges)
            )
        else:
            frames[frame_ix] = sm_compute_xymap_frame(
                ax_blocks[frame_ix], imsz, frame_ix, ctx, wcs, tranges
            )
    if pool is not None:
        pool.close()
        pool.join()
        frames = {ix: frame.get() for ix, frame in frames.items()}
    return write_backplane_movies(frames, imsz, ctx, tranges, wcs)


def make_backplanes(
    eclipse,
    band: GalexBand = "NUV",
    depth=None,
    leg=0,
    threads=None,
    burst=False,
    local='test_data',
    kind: Literal["xy", "dose"] = "xy",
    radius: int = 400,
    write: Optional[Mapping] = None,
    stop_after: Optional[str] = None,
    inline: bool = True,
    threshold: float = 0.75,
    star_size: float = 2,
    snippet: Optional[tuple] = None
):
    # noinspection PyTypeChecker
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
        stop_after=stop_after,
        snippet=snippet
    )
    # TODO: consider propagating this little hack upstream
    ctx.hdu_constructor_kwargs = dict(ctx.hdu_constructor_kwargs)
    ctx.inline = inline
    ctx.star_size = star_size
    ctx.threshold = threshold
    warn_bad_inline(ctx)
    if kind == "xy":
        return make_xymaps(ctx)
    elif kind == "dose":
        return make_dosemap(ctx, radius)
    raise ValueError("unrecognized backplane type")


def warn_bad_inline(ctx):
    combo = (
        ctx.write.get('array') is True, ctx.inline is True, ctx.burst is False
    )
    if all(combo):
        print("warning: inline array-writing only possible in burst mode.")


def check_inline_write(ctx):
    combo = (
        ctx.write.get('array') is True, ctx.inline is True, ctx.burst is True
    )
    if all(combo):
        return True
    return False


def write_or_return_arrays(maps, frame_ix, ctx, start_time, time_stamp, wcs=None, tranges=None):
    if check_inline_write(ctx):
        for k, v in maps.items():
            write_backplane_file(v, ctx, start_time, tranges, wcs, k, frame_ix, time_stamp)
        return
    elif ctx.write.get('array') is False:
        return
    return {
        name: pd.arrays.SparseArray(array.ravel())
        for name, array in maps.items()
    }


def dosemaps_just_for_timestamps(ctx: PipeContext, radius: int = 400):
    """ a function to make the list of time_stamps / frames generated by backplanes
    without actually using alll functionalities of backplanes. a little hacky.
     returns a dataframe of the time stamp and the modded timestamp for file naming"""
    components, start_time = load_for_dosemap(ctx()['photonfile'], radius, ctx.snippet)
    print('slicing position data into shared memory')
    #there's def a more efficient way to do this than using all that memory
    dose_blocks, tranges = send_to_shared_memory(components, ctx.depth)
    del components
    time_stamps = {}
    for frame_ix, trange in enumerate(tranges):
        # frame value is trange[0] - start time of eclipse
        time_stamps[frame_ix] = (trange[0], trange[0] - start_time.as_py())
    time_df = pd.DataFrame.from_dict(time_stamps, columns=['time', 'mod_time'], orient='index')
    return time_df


