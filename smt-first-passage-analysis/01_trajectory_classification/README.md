# 1 · Trajectory classification

`trajectory_classification.py` takes the raw TrackMate spots and tracks of each field of view (FOV), computes
features for every trajectory and classifies it. It is `260120-UTA-correctedcode.py`, unchanged.

## Input

A folder with matching pairs `<name>-spots.csv` / `<name>-tracks.csv` (or `_spots` / `_tracks`), **in pixels and
frames**. Rows that are not numbers, such as TrackMate 7's name, short-name and units rows, are removed when the
files are loaded.

## Settings (asked when the script starts)

| Prompt | Default | Used for the paper |
|---|---|---|
| Input folder | `data/input` | one condition folder |
| Output folder | `data/output` | |
| Frame interval (s) | 0.02 | 0.02 for MoeL and MoeLActin, 0.03 for MoeS and MoeSActin |
| Pixel to micron conversion | 0.097 | 0.097 |
| Minimum localizations for features | 5 | 5 |
| Minimum track duration for classification (s) | 0.2 | 0.18* |

\* Inferred from the saved results. In every FOV, the shortest classified track lasts 0.18 s and the longest
"Short" track with at least 5 localizations lasts 0.16 s (0.15 s at 30 ms), so the setting was between 0.16 and
0.18 s.

## Method

1. **Features for each track:**
   - minimum bounding circle radius (MBCR): the smallest maximum distance from the centroid or from any
     localization
   - distance between the MBC centre and the centre of mass (MBCC–CoM), and its ratio to MBCR
   - path length, tortuosity and elongation
   - MSD, and the diffusion coefficient from a linear fit to the first 4 lags, with or without intercept, whichever
     has the better R²
2. **Thresholds:** a 2-component Gaussian mixture model is fitted to the MBCR distribution, and then to the
   MBCC–CoM ratio of the high-MBCR tracks. Each threshold is the valley between the two components. If the GMM
   value is out of range (MBCR outside 0.05–0.3 µm, ratio outside 0.2–0.8), the literature values are used instead:
   118 nm and 0.35.
3. **Classification:**
   - *Short*: fewer localizations than the minimum, or shorter than the minimum duration
   - *Invalid*: MBCR below 0.001 µm, infinite tortuosity (the track ends where it started), or a NaN feature
   - *Immobile*: MBCR below the threshold
   - *Mobile*: MBCR above the threshold and MBCC–CoM ratio below its threshold
   - *Hybrid*: MBCR above the threshold and MBCC–CoM ratio above its threshold

## Output (for each FOV)

- `<name>_comprehensive_results.csv`: all features and the classification. This is the file the filter in step 2
  reads.
- `<name>_summary.csv`: key features
- `<name>_results.xlsx`: Summary, Statistics and All_Features sheets
- `<name>_analysis.png`: overview figure
- `<name>_report.txt`: thresholds and classification counts

Across all FOVs: `master_comprehensive_results.csv`, `master_analysis.png`, `master_distribution_plot.png` and
`master_report.txt`.
