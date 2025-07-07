import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.spatial.distance import cdist
import warnings
warnings.filterwarnings('ignore')
from pyarrow import ArrowInvalid
from tqdm import tqdm
from sklearn.cluster import DBSCAN
from  astropy.wcs import WCS

def indexify(cat,band='NUV'):
    return cat['ECLIPSE'].astype(str).str.zfill(5) + '_' + cat['LEG'].astype(str).str.zfill(2) + '_' + cat.index.astype(str).str.zfill(5) + '_' + band

def compute_separation_matrix(ra1, dec1, ra2, dec2):
    """
    Compute angular separations between two sets of coordinates efficiently.
    Uses astropy's SkyCoord separation matrix, which is vectorized and much faster.
    Returns a (len(ra1), len(ra2)) array of separations in arcseconds.
    """
    coords1 = SkyCoord(ra=ra1*u.degree, dec=dec1*u.degree)
    coords2 = SkyCoord(ra=ra2*u.degree, dec=dec2*u.degree)
    # Use astropy's separation matrix (broadcasts efficiently)
    # This returns a Quantity array of shape (len(coords1), len(coords2))
    sep_matrix = coords1[:, None].separation(coords2[None, :]).arcsec
    return sep_matrix

def imgpath_from_obsdata(eclipse, leg, band='nd', rootpath='data'):
    """
    Generate image file path from eclipse and leg information.
    
    Parameters:
    -----------
    eclipse : int
        Eclipse number
    leg : int  
        Leg number
    band : str
        Band identifier ('nd' for NUV, 'fd' for FUV)
    rootpath : str
        Root path to data directory
        
    Returns:
    --------
    str : Path to image file
    """
    estr = f'e{str(eclipse).zfill(5)}'
    return f'{rootpath}/{estr}/{estr}-{band}-ffull-b{str(leg).zfill(2)}-image-r.fits'


def make_wcs_from_header(header):
    """
    Create WCS object from FITS header.
    
    Parameters:
    -----------
    header : astropy.io.fits.Header
        FITS header containing WCS information
        
    Returns:
    --------
    astropy.wcs.WCS : WCS object
    """
    wcs = WCS(naxis=header['NAXIS'])
    wcs.wcs.cdelt = [header['CDELT1'], header['CDELT2']]
    wcs.wcs.ctype = [header['CTYPE1'], header['CTYPE2']]
    wcs.wcs.crpix = [header['CRPIX1'], header['CRPIX2']]
    wcs.wcs.crval = [header['CRVAL1'], header['CRVAL2']]
    return wcs


def get_catalog_sources_in_image(catalog, eclipse, leg, max_separation_arcmin=30.0):
    """
    Get all catalog sources that appear in the same eclipse/leg observation.
    
    Parameters:
    -----------
    catalog : pandas.DataFrame
        Source catalog with RA, DEC, ECLIPSE, LEG columns
    eclipse : int
        Eclipse number
    leg : int
        Leg number  
    max_separation_arcmin : float
        Maximum separation to include sources (arcminutes)
        
    Returns:
    --------
    pandas.DataFrame : Filtered catalog for this observation
    """
    # First filter by eclipse and leg
    obs_sources = catalog[
        (catalog['ECLIPSE'] == eclipse) & 
        (catalog['LEG'] == leg)
    ].copy()
    
    return obs_sources

def plot_source_thumbnail(catalog_row, catalog=None, band='NUV', rootpath='data', 
                         thumbnail_size_arcsec=300, figsize=(16, 8),
                         show_apertures=True, aperture_radii=[9.0, 17.5],
                         save_path=None, show_plot=True):
    eclipse = int(catalog_row['ECLIPSE'])
    leg = int(catalog_row['LEG'])
    source_ra = float(catalog_row['RA'])
    source_dec = float(catalog_row['DEC'])
    source_id = catalog_row['GLCAT_VISIT_ID']

    img_path = imgpath_from_obsdata(eclipse, leg, band[0].lower()+'d', rootpath)
    img_data = pdr.read(img_path)
    img_data.load('CNT')
    image = img_data['CNT']
    header = img_data['CNT_HEADER']
    
    wcs = make_wcs_from_header(header)
    source_coord = SkyCoord(ra=source_ra*u.degree, dec=source_dec*u.degree, frame='icrs')
    source_pixel = wcs.wcs_world2pix([[source_ra, source_dec]], 1)[0]
    source_x, source_y = source_pixel

    # get the x,y position of all stars in the catalog
    catalog_positions = wcs.wcs_world2pix(
            list(zip(catalog['RA'], catalog['DEC'])), 1
        )

    visit_positions = wcs.wcs_world2pix(
            list(zip(catalog[(catalog['ECLIPSE']==eclipse) & (catalog['LEG']==leg)]['RA'].values,
                    catalog[(catalog['ECLIPSE']==eclipse) & (catalog['LEG']==leg)]['DEC'].values)),
                    1
    )

    # Calculate thumbnail boundaries
    thumbnail_size_pixels = thumbnail_size_arcsec / 1.5  # Assuming ~1.5"/pixel
    x_min = max(0, int(source_y - thumbnail_size_pixels))
    x_max = min(image.shape[0], int(source_y + thumbnail_size_pixels))
    y_min = max(0, int(source_x - thumbnail_size_pixels))
    y_max = min(image.shape[1], int(source_x + thumbnail_size_pixels))

    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1], figure=fig)
    
    # Full frame image (left panel)
    ax1 = fig.add_subplot(gs[0], projection=wcs.celestial)

    ax1.imshow(ZScaleInterval()(image), origin='lower', cmap='Greys_r')
    thumbnail_rect = Rectangle(
        (y_min, x_min), y_max - y_min, x_max - x_min,
        linewidth=1, edgecolor='yellow', facecolor='none', linestyle='-'
    )
    ax1.add_patch(thumbnail_rect)

    radius_pixels = 30 / 1.5  # Convert to pixels
    circle = Circle((source_x, source_y), radius_pixels,
                linewidth=1, edgecolor='y', facecolor='none', 
                linestyle='--', alpha=0.8)
    ax1.add_patch(circle)

    if len(visit_positions) > 0:
        ax1.scatter(visit_positions[:, 0], visit_positions[:, 1], 
                c='b', s=4, marker='.', alpha=0.9, edgecolors='none')

    ax1.coords[0].set_axislabel('RA')
    ax1.coords[1].set_axislabel('Dec')
    ax1.coords[0].set_ticks(number=5)
    ax1.coords[1].set_ticks(number=5)
    ax1.coords[0].tick_params(labelsize=8)
    ax1.coords[1].tick_params(labelsize=8)

    ax1.set_title(f'{band} Full Frame', 
                fontsize=10, pad=10)

    # Thumbnail image (right panel)
    ax2 = fig.add_subplot(gs[1], projection=wcs.celestial)
    thumbnail = image[x_min:x_max, y_min:y_max]
    ax2.imshow(ZScaleInterval()(thumbnail), origin='lower', cmap='Greys_r',
                extent=(y_min, y_max, x_min, x_max))
    
    # Add yellow border around the thumbnail panel to match the box in full frame
    thumbnail_rect = Rectangle(
        (y_min, x_min), y_max - y_min, x_max - x_min,
        linewidth=4, edgecolor='yellow', facecolor='none', linestyle='-'
    )
    ax2.add_patch(thumbnail_rect)

    for i, radius in enumerate(aperture_radii):
        radius_pixels = radius / 1.5
        linestyle = '-' if i == 0 else '--'
        linewidth = 1# if i == 0 else 1
        circle = Circle((source_x, source_y), radius_pixels,
                    linewidth=linewidth, edgecolor='y', facecolor='none', 
                    linestyle=linestyle, alpha=0.9)
        ax2.add_patch(circle)

    thumb_mask = (
        (catalog_positions[:, 0] >= y_min) & (catalog_positions[:, 0] <= y_max) &
        (catalog_positions[:, 1] >= x_min) & (catalog_positions[:, 1] <= x_max)
    ) if len(catalog_positions) > 0 else []

    if len(catalog_positions) > 0 and np.any(thumb_mask):
        thumb_sources = catalog_positions[thumb_mask]
        ax2.scatter(thumb_sources[:, 0], thumb_sources[:, 1], 
                    c='yellow', s=12, marker='x', alpha=0.9)#, edgecolors='none')


    # ax2.scatter(catalog['RA'],catalog['DEC'],marker='x',color='y',extent=(y_min, y_max, x_min, x_max))

    ax2.coords[0].set_axislabel('RA') 
    ax2.coords[1].set_axislabel('Dec')
    ax2.coords[0].set_ticks(number=4)
    ax2.coords[1].set_ticks(number=4)
    ax2.coords[0].tick_params(labelsize=8)
    ax2.coords[1].tick_params(labelsize=8)

    ax2.set_title(f'Source Detail', fontsize=10, pad=10)
    
    # Add overall figure title
    fig.suptitle(f'{source_id}', fontsize=10, y=0.9)

    plt.show()
