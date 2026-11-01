# =============================================================================
# Actin Cortex Pore Size Analysis
# Based on: Eghiaian et al., Biophysical Journal 108(6), 2015
#
# Pipeline (Steps 3–5):
#   3a. Bandpass filtering
#   3b. Watershed segmentation
#   4.  Binary skeletonization
#   5.  Domain / pore size calculation
#
# Input:  folder of FOV*.tif AFM height images
# Output: filtered TIFs, label TIFs, skeleton TIFs, CSVs, histogram
# =============================================================================

# ── Install dependencies (run once in terminal) ──────────────────────────────
# pip install scikit-image tifffile numpy scipy matplotlib pandas

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import tifffile
import warnings
warnings.filterwarnings('ignore')

from scipy import ndimage
from scipy.ndimage import label as scipy_label
from scipy.stats import norm
from skimage import filters, morphology
from skimage.filters import threshold_otsu
from skimage.morphology import skeletonize
from skimage.measure import regionprops_table
from skimage.segmentation import watershed
from skimage.feature import peak_local_max


# =============================================================================
# USER PROMPTS — folder paths
# =============================================================================

while True:
    INPUT_FOLDER = input("Enter path to input folder containing FOV*.tif files: ").strip()
    if os.path.isdir(INPUT_FOLDER):
        break
    print(f"  Folder not found: '{INPUT_FOLDER}' — please try again.")

while True:
    OUTPUT_FOLDER = input("Enter path to output folder (will be created if needed): ").strip()
    try:
        os.makedirs(OUTPUT_FOLDER, exist_ok=True)
        break
    except Exception as e:
        print(f"  Could not create folder: {e} — please try again.")


# =============================================================================
# SETTINGS — edit these if needed
# =============================================================================

PIXEL_SIZE_NM = 107.50     # AFM pixel size in nanometres (15 nm in the paper)

# Bandpass filter cutoffs (in nm)
LOW_SIGMA_NM  = 10.0     # removes noise smaller than this
HIGH_SIGMA_NM = 500.0    # removes background variation larger than this

# Watershed settings
WATERSHED_COMPACTNESS = 0.001   # 0 = classic; small value improves compactness
MIN_DOMAIN_SIZE_NM2   = 500     # discard domains smaller than this (nm²)

SHOW_PLOTS = True        # False for headless / batch runs
COLORMAP   = 'afmhot'   # colormap for AFM height images


# =============================================================================
# DERIVED CONSTANTS (do not edit)
# =============================================================================

LOW_SIGMA_PX  = LOW_SIGMA_NM  / PIXEL_SIZE_NM
HIGH_SIGMA_PX = HIGH_SIGMA_NM / PIXEL_SIZE_NM
MIN_DOMAIN_PX = int(MIN_DOMAIN_SIZE_NM2 / (PIXEL_SIZE_NM ** 2))

folders = {
    'filtered'  : os.path.join(OUTPUT_FOLDER, '1_filtered'),
    'segmented' : os.path.join(OUTPUT_FOLDER, '2_segmented'),
    'skeleton'  : os.path.join(OUTPUT_FOLDER, '3_skeleton'),
    'overlay'   : os.path.join(OUTPUT_FOLDER, '4_overlay'),
    'results'   : os.path.join(OUTPUT_FOLDER, '5_results'),
}
for p in folders.values():
    os.makedirs(p, exist_ok=True)

# Collect FOV*.tif files, sorted by the numeric suffix
tif_files = sorted(
    glob.glob(os.path.join(INPUT_FOLDER, 'FOV*.tif')) +
    glob.glob(os.path.join(INPUT_FOLDER, 'FOV*.TIF')),
    key=lambda f: int(''.join(filter(str.isdigit,
                                     os.path.splitext(os.path.basename(f))[0]
                                     .replace('FOV', ''))) or '0')
)

if not tif_files:
    raise FileNotFoundError(
        f"No FOV*.tif files found in: {INPUT_FOLDER}\n"
        "Check that files are named FOV1.tif, FOV2.tif, etc."
    )

print(f"\nFound {len(tif_files)} FOV*.tif file(s):")
for f in tif_files:
    print(f"  {os.path.basename(f)}")
print(f"\nPixel size : {PIXEL_SIZE_NM} nm")
print(f"Bandpass   : {LOW_SIGMA_NM}–{HIGH_SIGMA_NM} nm  "
      f"({LOW_SIGMA_PX:.2f}–{HIGH_SIGMA_PX:.2f} px)")
print(f"Min domain : {MIN_DOMAIN_SIZE_NM2} nm²  ({MIN_DOMAIN_PX} px)")
print(f"Output     : {OUTPUT_FOLDER}\n")


# =============================================================================
# FUNCTIONS
# =============================================================================

def load_tif(path):
    """Load a TIFF file and return a 2-D float64 array."""
    img = tifffile.imread(path)
    if img.ndim == 3:
        img = img[0]
    elif img.ndim > 3:
        img = img[0, 0]
    return img.astype(np.float64)


def bandpass_filter(img, low_sigma, high_sigma):
    """
    Step 3a — Bandpass filter.
    Replicates the ImageJ bandpass: Gaussian(low) - Gaussian(high).
    Returns float64 array normalised to [0, 1].
    """
    smooth_low  = filters.gaussian(img, sigma=low_sigma,  preserve_range=True)
    smooth_high = filters.gaussian(img, sigma=high_sigma, preserve_range=True)
    bp = smooth_low - smooth_high
    bp -= bp.min()
    if bp.max() > 0:
        bp /= bp.max()
    return bp


def watershed_segment(filtered_img, compactness=0.001, min_size_px=5):
    """
    Step 3b — Grayscale watershed segmentation.

    1. Otsu threshold  →  binary mask (filaments = True)
    2. Distance transform on pore (dark) regions
    3. Local maxima inside pores as seeds
    4. Watershed on inverted distance map
    5. Remove small domains

    Returns
    -------
    labels      : int32 label image (0 = filament / background)
    binary_mask : bool array, True where filaments are detected
    """
    thresh      = threshold_otsu(filtered_img)
    binary_mask = filtered_img >= thresh
    pore_mask   = ~binary_mask
    distance    = ndimage.distance_transform_edt(pore_mask)

    coords = peak_local_max(
        distance,
        min_distance=max(2, int(LOW_SIGMA_PX)),
        labels=pore_mask.astype(int)
    )
    seed_mask = np.zeros(distance.shape, dtype=bool)
    if coords.size > 0:
        seed_mask[coords[:, 0], coords[:, 1]] = True
    markers, _ = scipy_label(seed_mask)

    labels = watershed(
        -distance,
        markers=markers,
        mask=pore_mask,
        compactness=compactness
    )
    labels = morphology.remove_small_objects(labels, min_size=min_size_px)
    return labels.astype(np.int32), binary_mask


def build_skeleton(binary_mask):
    """Step 4 — Thin the filament mask to 1-pixel-wide lines."""
    return skeletonize(binary_mask)


def measure_domains(labels, pixel_size_nm):
    """
    Step 5 — Pore size calculation.
    pore_size_nm = sqrt(area_px * pixel_size_nm^2)
    Returns a DataFrame with one row per domain.
    """
    props = regionprops_table(
        labels,
        properties=('label', 'area', 'centroid',
                    'axis_major_length', 'axis_minor_length',
                    'eccentricity', 'perimeter')
    )
    df = pd.DataFrame(props)
    df.rename(columns={
        'centroid-0': 'centroid_row_px',
        'centroid-1': 'centroid_col_px'
    }, inplace=True)

    px_area = pixel_size_nm ** 2
    df['area_px']       = df['area']
    df['area_nm2']      = df['area'] * px_area
    df['pore_size_nm']  = np.sqrt(df['area_nm2'])   # key formula
    df['major_axis_nm'] = df['axis_major_length'] * pixel_size_nm
    df['minor_axis_nm'] = df['axis_minor_length'] * pixel_size_nm
    df['aspect_ratio']  = np.where(
        df['axis_minor_length'] > 0,
        df['axis_major_length'] / df['axis_minor_length'],
        np.nan
    )
    df.drop(columns=['area', 'axis_major_length', 'axis_minor_length'],
            inplace=True)
    return df


# =============================================================================
# BATCH PROCESSING
# =============================================================================

all_domains  = []
summary_rows = []

for tif_path in tif_files:
    base = os.path.splitext(os.path.basename(tif_path))[0]
    print(f"Processing: {base}.tif")

    # Load
    raw = load_tif(tif_path)
    print(f"  Size: {raw.shape[1]} x {raw.shape[0]} px")

    # Step 3a — bandpass filter
    filtered        = bandpass_filter(raw, LOW_SIGMA_PX, HIGH_SIGMA_PX)
    filtered_uint16 = (filtered * 65535).astype(np.uint16)
    tifffile.imwrite(
        os.path.join(folders['filtered'], f"{base}_filtered.tif"),
        filtered_uint16, photometric='minisblack'
    )
    print(f"  Saved: 1_filtered/{base}_filtered.tif")

    # Step 3b — watershed segmentation
    labels, binary_mask = watershed_segment(
        filtered, WATERSHED_COMPACTNESS, MIN_DOMAIN_PX
    )
    n_domains = labels.max()
    print(f"  Domains found: {n_domains}")

    tifffile.imwrite(
        os.path.join(folders['segmented'], f"{base}_labels.tif"),
        labels.astype(np.int32), photometric='minisblack'
    )
    tifffile.imwrite(
        os.path.join(folders['segmented'], f"{base}_binary.tif"),
        binary_mask.astype(np.uint8) * 255, photometric='minisblack'
    )
    print(f"  Saved: 2_segmented/{base}_labels.tif  _binary.tif")

    # Step 4 — skeletonize
    skeleton   = build_skeleton(binary_mask)
    skel_uint8 = skeleton.astype(np.uint8) * 255
    tifffile.imwrite(
        os.path.join(folders['skeleton'], f"{base}_skeleton.tif"),
        skel_uint8, photometric='minisblack'
    )
    print(f"  Saved: 3_skeleton/{base}_skeleton.tif")

    # Step 5 — measure domains
    if n_domains == 0:
        print("  WARNING: no domains detected — skipping.\n")
        continue

    df = measure_domains(labels, PIXEL_SIZE_NM)
    df.insert(0, 'image', base)

    csv_path = os.path.join(folders['results'], f"{base}_domains.csv")
    df.to_csv(csv_path, index=False, float_format='%.3f')
    print(f"  Saved: 5_results/{base}_domains.csv  ({len(df)} domains)")

    all_domains.append(df)
    fov_mode_diam = pd.cut(df['pore_size_nm'], bins=30).value_counts().idxmax().mid
    fov_mode_area = pd.cut(df['area_nm2'],     bins=30).value_counts().idxmax().mid
    summary_rows.append({
        'image'            : base,
        'n_domains'        : len(df),
        'pore_size_mean_nm': df['pore_size_nm'].mean(),
        'pore_size_std_nm' : df['pore_size_nm'].std(),
        'pore_size_med_nm' : df['pore_size_nm'].median(),
        'pore_size_mode_nm': fov_mode_diam,        # <-- add
        'pore_size_min_nm' : df['pore_size_nm'].min(),
        'pore_size_max_nm' : df['pore_size_nm'].max(),
        'area_mean_nm2'    : df['area_nm2'].mean(),
        'area_mode_nm2'    : fov_mode_area,         # <-- add
        'aspect_ratio_mean': df['aspect_ratio'].mean(),
    })

    # 4-panel QC figure
    fig, axes = plt.subplots(1, 4, figsize=(18, 4))
    fig.suptitle(base, fontsize=11)

    axes[0].imshow(raw,      cmap=COLORMAP); axes[0].set_title('Raw AFM')
    axes[1].imshow(filtered, cmap=COLORMAP); axes[1].set_title('Bandpass filtered')

    label_rgb = plt.cm.tab20b((labels % 20) / 20.0)[..., :3]
    label_rgb[labels == 0] = 0
    label_rgb[skeleton]    = [1, 1, 1]
    axes[2].imshow(label_rgb); axes[2].set_title('Domains + skeleton')

    axes[3].hist(df['pore_size_nm'], bins=30,
                 color='steelblue', edgecolor='white', linewidth=0.4)
    axes[3].axvline(df['pore_size_nm'].mean(),   color='red',    lw=1.5, label='mean')
    axes[3].axvline(df['pore_size_nm'].median(), color='orange', lw=1.5,
                    ls='--', label='median')
    axes[3].set_xlabel('Pore size (nm)')
    axes[3].set_ylabel('Count')
    axes[3].set_title('Pore size distribution')
    axes[3].legend(fontsize=8)
    for ax in axes[:3]:
        ax.axis('off')
    plt.tight_layout()

    fig_path = os.path.join(folders['overlay'], f"{base}_overview.png")
    plt.savefig(fig_path, dpi=150, bbox_inches='tight')
    if SHOW_PLOTS:
        plt.show()
    plt.close()
    print(f"  Saved: 4_overlay/{base}_overview.png\n")


# =============================================================================
# AGGREGATE RESULTS
# =============================================================================

if not all_domains:
    print("No domains measured — check input folder and image content.")
else:
    df_summary   = pd.DataFrame(summary_rows)
    summary_path = os.path.join(folders['results'], 'summary_per_image.csv')
    df_summary.to_csv(summary_path, index=False, float_format='%.3f')
    print(f"Saved: {summary_path}")
    print(df_summary.round(2).to_string(index=False))

    df_all   = pd.concat(all_domains, ignore_index=True)
    all_path = os.path.join(folders['results'], 'all_domains_pooled.csv')
    df_all.to_csv(all_path, index=False, float_format='%.3f')
    print(f"\nSaved pooled data: {all_path}")
    print(f"Total domains: {len(df_all)} across {len(all_domains)} image(s)")

    ps = df_all['pore_size_nm']
    global_mode_diam = pd.cut(ps, bins=40).value_counts().idxmax().mid
    global_mode_area = pd.cut(df_all['area_nm2'], bins=40).value_counts().idxmax().mid
    print(f"\n── Global pore size statistics ──")
    print(f"  Mean  +/- SD : {ps.mean():.1f} +/- {ps.std():.1f} nm")
    print(f"  Median       : {ps.median():.1f} nm")
    print(f"  Mode (diam)  : {global_mode_diam:.1f} nm")
    print(f"  Mode (area)  : {global_mode_area:.1f} nm²")
    print(f"  Range        : {ps.min():.1f} – {ps.max():.1f} nm")
    print(f"  (Paper: 96.5 +/- 72.2 nm, 6070 domains, 124 frames)")

    # Global histogram with Gaussian fit (replicates Fig. 2E)
    fig, ax = plt.subplots(figsize=(7, 4))
    n, bins, _ = ax.hist(
        ps, bins=40, range=(0, 400),
        color='#888888', edgecolor='white', linewidth=0.5,
        label=f'n = {len(ps):,} domains'
    )
    mu, sigma = ps.mean(), ps.std()
    x         = np.linspace(0, 400, 300)
    bin_width = bins[1] - bins[0]
    gauss     = norm.pdf(x, mu, sigma) * len(ps) * bin_width
    ax.plot(x, gauss, 'r-', lw=2,
            label=f'Gaussian fit\nmu={mu:.1f} nm, sigma={sigma:.1f} nm')
    ax.axvline(mu, color='red', lw=1, ls='--', alpha=0.7)
    ax.axvline(global_mode_diam, color='green', lw=1.5, ls=':',
               label=f'Mode = {global_mode_diam:.1f} nm\n(area mode = {global_mode_area:.0f} nm²)')
    ax.set_xlabel('Domain (pore) size (nm)', fontsize=12)
    ax.set_ylabel('Events', fontsize=12)
    ax.set_title('Actin cortex pore size distribution', fontsize=12)
    ax.legend(fontsize=9)
    ax.set_xlim(0, 400)
    plt.tight_layout()

    hist_path = os.path.join(folders['results'], 'pore_size_histogram.png')
    plt.savefig(hist_path, dpi=150, bbox_inches='tight')
    print(f"\nSaved: {hist_path}")
    if SHOW_PLOTS:
        plt.show()
    plt.close()

print("\nDone.")
