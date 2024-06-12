import pandas as pd
import numpy as np
import fast_histogram as fh
import astropy.wcs
import os
from gPhoton.coadd import make_shared_wcs, project_to_shared_wcs, bin_projected_weights
from astropy.visualization import ZScaleInterval
import matplotlib.pyplot as plt
import imageio

def project_pixels_to_new_wcs(x, y, from_wcs, to_wcs):
    ra, dec = from_wcs.pixel_to_world_values(x, y)
    return to_wcs.wcs_world2pix(ra, dec, 1)

def imsz_from_wcs(reference_wcs):
    return (
        int((reference_wcs.wcs.crpix[1] - 0.5) * 2),
        int((reference_wcs.wcs.crpix[0] - 0.5) * 2),
    )

def bin_projected_weights(x, y, weights, imsz):
    binned = fh.histogram2d(
        y - 0.5,
        x - 0.5,
        bins=imsz,
        range=([[0, imsz[0]], [0, imsz[1]]]),
        weights=weights,
    )
    return binned

def wcs_imsz(system: astropy.wcs.WCS):
    """
    image size associated with a WCS object. WARNING: not universally
    applicable! works if and only if the reference pixel is at the center of
    the image.
    """
    return (
        int((system.wcs.crpix[1] - 0.5) * 2),
        int((system.wcs.crpix[0] - 0.5) * 2),
    )


def project_to_shared_wcs(
    image: np.ndarray,
    from_wcs: astropy.wcs.WCS,
    to_wcs: astropy.wcs.WCS,
    ):
    y_ix, x_ix = np.nonzero(image)
    ra_input, dec_input = from_wcs.pixel_to_world_values(x_ix, y_ix)
    x_shared, y_shared = to_wcs.wcs_world2pix(ra_input, dec_input, 1)
    return {
        "x": x_shared,
        "y": y_shared,
        "weight": image[y_ix, x_ix],
        "imsz":wcs_imsz(to_wcs),
    }


def project_image_to_new_wcs(image: np.ndarray,
                             from_wcs, to_wcs):
    proj = project_to_shared_wcs(image,from_wcs,to_wcs)#,1,apply_mask=apply_mask)
    return bin_projected_weights(proj['x'],proj['y'],proj['weight'],proj['imsz'])

def make_band_gif(data,fn_root,temp_dir,output_dir,
                  overwrite = True,cmap="Greys_r",fps=1,apply_mask=True,
                  catalog=None,catalog_wcs=None,figsize=(10,10)):
    gif_fn = f'{output_dir}{fn_root}.gif'
    if not (catalog is None) == (catalog_wcs is None):
        raise ValueError('catalog and catalog_wcs must both best specified or neither')
    
    if os.path.exists(gif_fn) and not overwrite:
        return gif_fn
    
    shared_wcs = make_shared_wcs((data['NUV']['wcs'],data['FUV']['wcs']))

    for band in ['NUV','FUV']:
        frame_fn = f'{temp_dir}{fn_root}-{band}-catqa.png'
        if os.path.exists(frame_fn) and not overwrite:
            continue

        pimg = project_image_to_new_wcs(data[band]['count'],
                                        data[band]['wcs'],
                                        shared_wcs)
        
        plt.figure(figsize=figsize);
        plt.xticks([])
        plt.yticks([])
        plt.tight_layout()
        plt.imshow(ZScaleInterval()(pimg),cmap=cmap,origin="lower")

        if not catalog is None:
            x_shared,y_shared = project_pixels_to_new_wcs(catalog['xcentroid'].values,
                                                          catalog['ycentroid'].values,
                                                          catalog_wcs,shared_wcs)
            plt.plot(x_shared,y_shared,'y.',markersize=2)
            
        plt.savefig(frame_fn)
        plt.close()
    if os.path.exists(f'{output_dir}{fn_root}.gif') and not overwrite:
        pass
    else:
        with imageio.get_writer(gif_fn, mode='I', fps=fps, loop=0) as writer:
            for band in ['NUV','FUV']:
                frame_fn = f'{temp_dir}{fn_root}-{band}-catqa.png'
                image = imageio.imread(f'{frame_fn}')
                writer.append_data(image)
    return gif_fn

def crossmatch_catalogs(master_table, match_table, match_radius=3*4.17e-4):
    match_catalog = SkyCoord(ra=match_table['ra'].values*u.degree,
                             dec=match_table['dec'].values*u.degree)
    master_catalog = SkyCoord(ra=master_table['ra'].values*u.degree,
                              dec=master_table['dec'].values*u.degree)
    catalog_ix, d2d, d3d = match_catalog.match_to_catalog_sky(master_catalog)
    mask = np.ones(np.array(d2d).shape)
    mask[np.where(np.array(d2d)>match_radius)] = 0
    d2d_masked = np.ma.array(d2d,mask=mask)
    catalog_ix_masked = np.ma.array(catalog_ix,mask=mask)
    return d2d_masked, catalog_ix_masked