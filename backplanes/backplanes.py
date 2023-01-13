from functools import partial
from itertools import product
from multiprocessing import Pool
from pathlib import Path
from typing import Optional

import astropy.io.fits
import fast_histogram as fh
import numpy as np
import pandas as pd
import sh
from more_itertools import windowed
from pyarrow import parquet
from scipy.stats import binned_statistic_2d

from gPhoton.moviemaker._steps import (
    select_on_detector,
    generate_wcs_components,
    slice_into_memory
)
from gPhoton.pretty import print_inline
from gPhoton.reference import eclipse_to_paths
from gPhoton.sharing import (
    reference_shared_memory_arrays
)
from gPhoton.types import GalexBand
from gPhoton.vorpal import between


def components_to_shared_memory(components, depth):
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


def load_for_dosemap(photonfile, radius=400):
    phot = parquet.read_table(
        photonfile, columns=['t', 'col', 'row', 'detrad']
    )
    phot = select_on_detector(phot, radius)
    return {
        't': phot['t'].to_numpy(),
        'x': phot['row'].to_numpy() * 4,
        'y': phot['col'].to_numpy() * 4
    }


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
    )


def sm_make_dosemap(block_info, radius, frame_ix):
    _, axes = reference_shared_memory_arrays(block_info)
    print(f'integrating frame {frame_ix} dosemap')
    return {'dose': dosemap_frame(axes, radius)}


def sm_compute_dosemap_frame(block_info, imsz, frame_ix):
    maps = sm_make_dosemap(block_info, imsz, frame_ix)
    for _, block in reference_shared_memory_arrays(
        block_info, fetch=False
    )[0].items():
        block.close()
        block.unlink()
    return {
        name: pd.arrays.SparseArray(array.ravel())
        for name, array in maps.items()
    }


def make_dosemap(
    eclipse: int,
    band: GalexBand = 'NUV',
    leg: int = 0,
    depth: Optional[int] = None,
    burst: bool = False,
    root: str = 'data',
    threads: Optional[int] = None,
    radius: int = 400
):
    paths = eclipse_to_paths(eclipse, band, leg=leg, root=root)
    components = load_for_dosemap(paths['photonfile'], radius)
    if depth is None:
        maps = dosemap_frame(components, radius)
        return write_backplane_image(maps, eclipse, band, leg, root)
    else:
        print('slicing position data into shared memory')
        dose_blocks, tranges = components_to_shared_memory(components, depth)
        del components
    frames = {}
    pool = Pool(threads) if threads is not None else None
    for frame_ix, trange in enumerate(tranges):
        if pool is not None:
            frames[frame_ix] = pool.apply_async(
                sm_make_dosemap, (dose_blocks[frame_ix], radius, frame_ix)
            )
        else:
            frames[frame_ix] = sm_compute_dosemap_frame(
                dose_blocks[frame_ix], radius, frame_ix
            )
    if pool is not None:
        pool.close()
        pool.join()
        frames = {ix: frame.get() for ix, frame in frames.items()}
    ranges = dosemap_ranges(radius)
    imsz = [ranges[0][1] - ranges[0][0]] * 2
    return write_backplane_movies(
        frames, imsz, eclipse, band, leg, depth, burst, root
    )


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
        range=([[0, imsz[0]], [0, imsz[1]]]).append,
        statistic=stat
    )[0].astype(np.float32)


def sm_make_xymaps(block_info, imsz, frame_ix):
    _, arrays = reference_shared_memory_arrays(block_info)
    xymaps = {}
    for ax, stat in product(('x', 'y'), ('mean', 'std')):
        print_inline(f'integrating {frame_ix} {ax}_{stat}')
        xymaps[f'{ax}_{stat}'] = sm_binner(arrays, ax, stat, imsz)
    return xymaps


def sm_compute_xymap_frame(block_info, imsz, frame_ix):
    maps = sm_make_xymaps(block_info, imsz, frame_ix)
    for _, block in reference_shared_memory_arrays(
        block_info, fetch=False
    )[0].items():
        block.close()
        block.unlink()
    return {
        name: pd.arrays.SparseArray(array.ravel())
        for name, array in maps.items()
    }


def write_backplane_file(
    image, root, eclipse, band, leg, depth, name, frame_ix='movie'
):
    basename = (
        f'e{eclipse}_{band[0].lower()}d_{leg}b_{depth}s_{frame_ix}_{name}.fits'
    )
    for ext in ('', '.gz'):
        if Path(root, basename + ext).exists():
            Path(basename).unlink()
    hdu = astropy.io.fits.PrimaryHDU(image)
    hdu.writeto(basename)
    sh.igzip(basename)
    Path(basename).unlink()


def write_backplane_image(maps, eclipse, band, leg, root):
    for name, image in maps.items():
        write_backplane_file(image, root, eclipse, band, leg, name, 'image')


def write_backplane_movies(
    frames, imsz, eclipse, band, leg, depth, burst, root
):
    title = (root, eclipse, band, leg, depth)
    if burst is False:
        maps = {k: [] for k in frames[0].keys()}
        for ix in list(frames.keys()):
            for name, array in frames[ix].items():
                maps[name].append(array)
            del frames[ix]
        for name in list(maps.keys()):
            print(f'writing {name}')
            write_backplane_file(
                np.dstack([f.to_dense().reshape(imsz) for f in maps[name]]),
                *title,
                name
            )
            del maps[name]
        return
    for ix in list(frames.keys()):
        for name, sparse in frames[ix].items():
            print(f'writing {name} (frame {ix})')
            write_backplane_file(
                sparse.to_dense().reshape(imsz),
                *title,
                name,
                str(ix).zfill(5)
            )
        del frames[ix]


def make_xymaps(
    eclipse: int,
    band: GalexBand = 'NUV',
    leg: int = 0,
    depth: Optional[int] = None,
    burst: bool = False,
    root: str = 'data',
    threads: Optional[int] = None
):
    paths = eclipse_to_paths(eclipse, band, leg=leg, root=root)
    print("loading photonlist and computing wcs")
    components, wcs, imsz = load_for_xymap(paths['photonfile'])
    if depth is None:
        maps = make_full_depth_xymaps(components, imsz)
        return write_backplane_image(maps, eclipse, band, leg, root)
    else:
        print('slicing position data into shared memory')
        ax_blocks, tranges = components_to_shared_memory(components, depth)
        del components
    frames = {}
    pool = Pool(threads) if threads is not None else None
    for frame_ix, trange in enumerate(tranges):
        if pool is not None:
            frames[frame_ix] = pool.apply_async(
                sm_compute_xymap_frame, (ax_blocks[frame_ix], imsz, frame_ix)
            )
        else:
            frames[frame_ix] = sm_compute_xymap_frame(
                ax_blocks[frame_ix], imsz, frame_ix
            )
    if pool is not None:
        pool.close()
        pool.join()
        frames = {ix: frame.get() for ix, frame in frames.items()}
    return write_backplane_movies(
        frames, imsz, eclipse, band, leg, depth, burst, root
    )
