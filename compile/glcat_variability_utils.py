import pandas as pd
from rich import print
import numpy as np
import matplotlib.pyplot as plt
#from quickbin import bin2d
from scipy.stats import binned_statistic_2d
import pdr
from astropy.visualization import ZScaleInterval
from astropy import wcs as pywcs
import warnings
from pyarrow import ArrowInvalid
from tqdm import tqdm
from scipy.stats import anderson
from sklearn.cluster import DBSCAN
from scipy.spatial import ConvexHull
from matplotlib.patches import Polygon

import pyarrow.parquet as pq
table = pq.read_table('data/e23456/e23456-nd-f0120-b00-movie-photom-9_0-base.parquet').to_pandas()
data = pdr.read('data/e23456/e23456-nd-f0120-b00-movie-r.fits')
data.load('CNT')
image = pdr.read('data/e23456/e23456-nd-ffull-b00-image-r.fits')

nbins = data['CNT_HEADER']['N_FRAME']
expt = np.array([data['CNT_HEADER'][f'EXPT_{n}'] for n in range(nbins)])
t = np.arange(nbins)

def extract_photometry_data(table, nbins, expt):
    """
    Extract aperture counts and flags from table and calculate CPS and errors.
    
    Parameters:
    -----------
    table : pandas.DataFrame
        Table containing aperture_sum_{i} and artifact_flag_{i} columns
    nbins : int
        Number of time bins
    expt : array-like
        Exposure times for each bin
        
    Returns:
    --------
    tuple
        (aper_cnts_all, flags_all, cps_all, cps_err_all, good_mask)
    """
    # Generate column names
    aperture_cols = [f'aperture_sum_{i}' for i in range(nbins)]
    flag_cols = [f'artifact_flag_{i}' for i in range(nbins)]

    # Extract data arrays
    aper_cnts_all = table[aperture_cols].values  # Shape: (n_sources, nbins)
    flags_all = table[flag_cols].values  # Shape: (n_sources, nbins)

    # Calculate CPS and errors for all sources at once
    cps_all = aper_cnts_all / expt[np.newaxis, :]  # Broadcast expt across sources
    cps_err_all = np.sqrt(aper_cnts_all) / expt[np.newaxis, :]

    # Create mask for good (unflagged) data
    good_mask = (flags_all == 0)
    
    return aper_cnts_all, flags_all, cps_all, cps_err_all, good_mask

def detect_variability_outliers(table, nbins, expt, sigma=3, min_outliers=1):
    """
    Detect variability outliers in photometric time series data using vectorized operations.
    
    Parameters:
    -----------
    table : pandas.DataFrame
        Table containing aperture_sum_{i} and artifact_flag_{i} columns for each time bin
    nbins : int
        Number of time bins
    expt : array-like
        Exposure times for each bin
    sigma : float, optional
        Sigma threshold for outlier detection (default: 3)
    min_outliers : int, optional
        Minimum number of outliers required to flag a source (default: 2)
        
    Returns:
    --------
    list
        List of [outlier_count, source_index] pairs for sources with sufficient outliers
    """
    # Extract and preprocess photometry data
    aper_cnts_all, flags_all, cps_all, cps_err_all, good_mask = extract_photometry_data(table, nbins, expt)

    # Count good measurements per source
    n_good_per_source = np.sum(good_mask, axis=1)
    valid_sources = n_good_per_source > 2

    # Create arrays with NaN for invalid measurements
    cps_masked = np.where(good_mask, cps_all, np.nan)
    cps_err_masked = np.where(good_mask, cps_err_all, np.nan)

    # Calculate upper and lower bounds for each measurement
    upper_bounds = cps_masked + cps_err_masked * sigma
    lower_bounds = cps_masked - cps_err_masked * sigma

    # Sort each row and get second min/max (ignoring NaNs)
    sorted_upper = np.sort(upper_bounds, axis=1)  # For second_min calculation
    sorted_lower = np.sort(lower_bounds, axis=1)  # For second_max calculation

    # Get second minimum and second maximum for each source
    # We need to account for NaNs which sort to the end
    second_min_indices = np.minimum(1, n_good_per_source - 1)  # Clamp to valid range
    second_max_indices = np.maximum(n_good_per_source - 2, 0)  # Clamp to valid range

    # Use advanced indexing to get the values
    row_indices = np.arange(len(table))
    second_mins = sorted_upper[row_indices, second_min_indices]
    second_maxs = sorted_lower[row_indices, second_max_indices]

    # Count outliers for all sources at once
    outliers_high = np.sum((lower_bounds > second_mins[:, np.newaxis]) & good_mask, axis=1)
    outliers_low = np.sum((upper_bounds < second_maxs[:, np.newaxis]) & good_mask, axis=1)

    # Total outlier counts
    outlier_counts = outliers_high + outliers_low

    # Zero out counts for invalid sources
    outlier_counts = np.where(valid_sources, outlier_counts, 0)

    # Filter sources with sufficient outliers
    outlier_mask = (valid_sources) & (outlier_counts >= min_outliers)
    outlier_indices = np.where(outlier_mask)[0]

    # Create outlier list in the same format as original
    outliers = [table.index[i] for i in outlier_indices]
    
    # second method using Anderson-Darling
    out = [anderson(cps[:-1]) for cps in cps_all]
    ix = np.where([o.statistic>o.critical_values.min() for o in out])

    vix = list(set(ix[0].tolist()+outliers))

    return vix

# Example usage:
# outlier = detect_variability_outliers(table, nbins, expt, sigma=3, min_outliers=2)
print(f"Found {len(outliers:=detect_variability_outliers(table, nbins, expt))} variable sources")
#plt.semilogy()
# plt.show()

# aper_cnts_all, flags_all, cps_all, cps_err_all, good_mask = extract_photometry_data(table, nbins, expt)

# out = [anderson(cps[:-1]) for cps in cps_all]
# ix = np.where([o.statistic>o.critical_values.min() for o in out])
# print(ix)

# vix = list(set(ix[0].tolist()+outliers))
# clustering = DBSCAN(eps=1/60,min_samples=1).fit(table.iloc[vix][['ra','dec']].values)
# clustering.labels_

# for i in set(clustering.labels_):
#     if sum(clustering.labels_==i)>1:
#         plt.figure()
#         plt.title(i)
#         plt.plot(table.iloc[np.array(vix)[clustering.labels_==i]]['ra'],
#                  table.iloc[np.array(vix)[clustering.labels_==i]]['dec'],'k.')

from catalog_visualization import make_wcs_from_header
from scipy.stats import anderson
from sklearn.cluster import DBSCAN
image.load('CNT',reload=True)
wcs = make_wcs_from_header(image['CNT_HEADER'])
plt.figure()
plt.imshow(ZScaleInterval()(image['CNT']), origin='lower', cmap='Greys_r')
plt.plot(table['xcenter'],table['ycenter'],'yx')
plt.plot(table.iloc[np.array(vix)]['xcenter'],
         table.iloc[np.array(vix)]['ycenter'],'bx')


###########################

# Run DBSCAN on all ra,dec positions in the table and plot convex hulls
# Run DBSCAN clustering on all sources in the table
print("Running DBSCAN clustering on all sources...")
# Use a reasonable epsilon for astronomical coordinates (e.g., 30 arcseconds = 30/3600 degrees)
eps_deg = 15/3600  # 30 arcseconds in degrees
min_samples = 3  # Minimum sources to form a cluster

clustering_all = DBSCAN(eps=eps_deg, min_samples=min_samples).fit(table[['ra','dec']].values)
n_clusters = len(set(clustering_all.labels_)) - (1 if -1 in clustering_all.labels_ else 0)
n_noise = list(clustering_all.labels_).count(-1)

print(f"Found {n_clusters} clusters and {n_noise} noise points")

# Create the plot
fig, ax = plt.subplots(figsize=(16, 12))

# Load and display the image
image.load('CNT', reload=True)
wcs = make_wcs_from_header(image['CNT_HEADER'])

# Display the image with proper scaling
img_data = image['CNT']
vmin, vmax = ZScaleInterval().get_limits(img_data)
ax.imshow(img_data, origin='lower', cmap='Greys_r', vmin=vmin, vmax=vmax)

# Convert ra,dec to pixel coordinates for plotting
source_coords = wcs.wcs_world2pix(table[['ra','dec']].values, 1)

# Define colors for clusters
colors = plt.cm.tab10(np.linspace(0, 1, min(10, n_clusters)))
if n_clusters > 10:
    # If more than 10 clusters, use a larger colormap
    colors = plt.cm.tab20(np.linspace(0, 1, min(20, n_clusters)))

# Plot each cluster and its convex hull
for cluster_id in set(clustering_all.labels_):
    if cluster_id == -1:  # Skip noise points
        continue
    
    # Get sources in this cluster
    cluster_mask = clustering_all.labels_ == cluster_id
    cluster_sources = table[cluster_mask]
    cluster_pixels = source_coords[cluster_mask]
    
    # does this cluster contain multiple variables?
    if not len(set(cluster_sources.index).intersection(outliers))>1:
        continue

    if len(cluster_sources) < 3:
        # Need at least 3 points for a convex hull
        continue
    
    # Get color for this cluster
    color_idx = cluster_id % len(colors)
    color = colors[color_idx]
    
    # Plot the sources in this cluster
    ax.scatter(cluster_pixels[:, 0], cluster_pixels[:, 1], 
              c=[color], s=20, alpha=0.8, 
              label=f'Cluster {cluster_id} ({len(cluster_sources)} sources)')
    
    # Calculate and plot convex hull
    try:
        hull = ConvexHull(cluster_pixels)
        hull_points = cluster_pixels[hull.vertices]
        
        # Create polygon for the convex hull
        hull_polygon = Polygon(hull_points, fill=False, 
                              edgecolor=color, linewidth=2, 
                              linestyle='--', alpha=0.7)
        ax.add_patch(hull_polygon)
        
        # Add cluster label at centroid
        centroid = np.mean(cluster_pixels, axis=0)
        ax.annotate(f'C{cluster_id}', centroid, 
                   xytext=(5, 5), textcoords='offset points',
                   fontsize=8, color=color, weight='bold',
                   bbox=dict(boxstyle='round,pad=0.2', 
                            facecolor='white', alpha=0.7))
        
    except Exception as e:
        continue
        # print(f"Could not create convex hull for cluster {cluster_id}: {e}")

# Plot noise points (if any)
# if n_noise > 0:
#     noise_mask = clustering_all.labels_ == -1
#     noise_pixels = source_coords[noise_mask]
#     ax.scatter(noise_pixels[:, 0], noise_pixels[:, 1], 
#               c='y', s=10, alpha=0.5, marker='x',
#               label=f'Noise ({n_noise} sources)')

# Set labels and title
ax.set_xlabel('X (pixels)', fontsize=12)
ax.set_ylabel('Y (pixels)', fontsize=12)
ax.set_title(f'DBSCAN Clustering of All Sources with Convex Hulls\n'
            f'Parameters: eps={eps_deg*3600:.1f}", min_samples={min_samples}\n'
            f'Found {n_clusters} clusters, {n_noise} noise points', 
            fontsize=14, pad=20)

# Add legend (limit to reasonable number of entries)
handles, labels = ax.get_legend_handles_labels()
if len(handles) > 15:
    # If too many clusters, just show summary
    ax.text(0.02, 0.98, f'{n_clusters} clusters found\nSee annotations for cluster IDs', 
           transform=ax.transAxes, fontsize=10, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
else:
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)

# Set aspect ratio and tight layout
ax.set_aspect('equal')
plt.tight_layout()
plt.show()

# Print cluster statistics
print(f"\nCluster Statistics:")
print(f"Total sources: {len(table)}")
print(f"Clustered sources: {len(table) - n_noise}")
print(f"Number of clusters: {n_clusters}")
if n_clusters > 0:
    cluster_sizes = [sum(clustering_all.labels_ == i) for i in set(clustering_all.labels_) if i != -1]
    print(f"Cluster size range: {min(cluster_sizes)} - {max(cluster_sizes)} sources")
    print(f"Average cluster size: {np.mean(cluster_sizes):.1f} sources")

