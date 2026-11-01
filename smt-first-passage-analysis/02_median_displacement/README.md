# 2 · Mobile track filter and median displacement

## `sort_mobile_tracks.py`: mobile track filter

This script selects the tracks used in steps 3–5. It is `Sortingmobiletracks.docx`, copied as text without changes.

### Input

A base folder with `MoeL/`, `MoeLActin/`, `MoeS/` and `MoeSActin/`. Each contains, for every FOV:

- `FOV<n>-spots.csv` and `FOV<n>-tracks.csv` in µm and seconds
- `FOV<n>_comprehensive_results.csv` from step 1. If it is missing, all tracks are used and a warning is printed.

### Settings

- **Asked when the script starts:** the base folder, and 1–3 duration thresholds in seconds. The manuscript FPT/FET
  figures used 1.0 s (`Rawdata-mobiletracks-1sec`). A 0.6 s set was also made (`Rawdata-mobiletracks-600msec`).
- **Set in the code:** the 0.2 µm displacement threshold, and the frame interval for each condition (0.02 s for
  MoeL and MoeLActin, 0.03 s for MoeS and MoeSActin).

### Method

For each FOV:

1. Keep the tracks labelled **Mobile** in the comprehensive results file.
2. Of those, keep the tracks whose **median frame-to-frame displacement is at least 0.2 µm**.
3. Of those, keep the tracks whose `TRACK_DURATION` is **at least each duration threshold**.

### Output

The output is written into each condition folder:

- `FOV<n>-0.2um-mobile spots.csv` / `… tracks.csv`: tracks after step 2 of the method
- `FOV<n>-0.2um-dur<t>s-mobile spots.csv` / `… tracks.csv`: tracks after step 3, for each threshold
- `<condition>-summary.csv`, and `all_conditions_summary.csv` in the base folder: track counts after each filter

Steps 3–5 look for files named `FOV<n>-spots.csv`. Copy the files for the threshold you want (e.g. `dur1s`) to
a new base folder with the same condition subfolders, and rename them to that pattern.

## `median_displacement.py`: displacement plots

This optional script plots the trajectories of one condition folder. It is `Median-displacement-gemini-code.docx`,
copied as text without changes.

### Input

A folder with `FOV*-spots.csv` files.

### Settings (asked when the script starts)

| Prompt | Default |
|---|---|
| Input directory | none |
| Output directory | none |
| Frame interval (s) | none. Asked for, but not used in the calculation |
| Max scale for spider plot axes (µm) | 2 |
| Minimum duration threshold (s) | 1.0 |

### Method

1. Keep tracks whose duration (last time minus first time) is at least the minimum.
2. Compute the frame-to-frame displacement and the median displacement of each track.
3. *Mobile* tracks are those with a median displacement of at least 0.2 µm (hard-coded).

### Output (for each FOV)

- `<FOV>_dist.png`: histogram and KDE of the median displacement per track, marking the mean, the median and the
  200 nm threshold
- `<FOV>_All_spatial.png` and `<FOV>_All_spider.png`: all tracks that pass the duration filter
- `<FOV>_Mobile_spatial.png`, `<FOV>_Mobile_spider.png` and `<FOV>_mobile_dynamics.png`: mobile tracks only
- `report.csv`: number of tracks per FOV that passed the duration filter, fell below 200 nm, or remained as mobile

The duration is computed from `POSITION_T`. Tracks lasting exactly the minimum can compute as slightly shorter
(e.g. 0.5999… s) and be dropped.
