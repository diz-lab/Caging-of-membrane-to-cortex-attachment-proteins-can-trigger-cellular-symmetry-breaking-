# Single-molecule tracking: first passage time analysis

Python pipeline for analysing single-molecule tracking (SMT) trajectories exported from
[TrackMate](https://imagej.net/plugins/trackmate/) (Fiji). It classifies and filters trajectories, then measures
first passage times (FPT), free exit times (FET), track lifetimes, end-to-end distances and path lengths. It
compares four conditions (MoeS, MoeL, MoeSActin, MoeLActin) and draws spider plots of the trajectories.

## Pipeline

```mermaid
flowchart LR
    A[TrackMate export<br/>pixels + frames] --> B[1 · Trajectory<br/>classification]
    B -->|FOV*_comprehensive_results.csv| C[2 · Mobile track filter<br/>+ median displacement]
    C -->|filtered FOV*-spots.csv| D[3 · First passage<br/>time]
    D --> E[4 · Bimodal<br/>FPT / FET fitting]
    E --> F[5 · Spider<br/>plots]
```

| Step | Folder | What it does |
|---|---|---|
| 1 | [`01_trajectory_classification`](01_trajectory_classification) | Computes per-track features (minimum bounding circle radius, MBCC–CoM ratio, MSD, diffusion coefficient, tortuosity, …) and classifies tracks as Immobile, Mobile, Hybrid, Short or Invalid. Thresholds are set with Gaussian-mixture models. |
| 2 | [`02_median_displacement`](02_median_displacement) | `sort_mobile_tracks.py` keeps the tracks labelled Mobile in step 1, then those with a median step of at least 200 nm, then those lasting at least the chosen duration (1 s for the manuscript FPT/FET figures), and saves the filtered spots and tracks. `median_displacement.py` plots the displacement distributions, trajectories and spider plots for each FOV. |
| 3 | [`03_first_passage_time`](03_first_passage_time) | First passage time at radii of 1, 2 and 3 µm, plus free exit time, lifetime, end-to-end distance and path length for each condition, with Gaussian fits. |
| 4 | [`04_bimodal_fpt_fet`](04_bimodal_fpt_fet) | The manuscript FPT/FET figures: more radii and thresholds, automatic choice between Gaussian, exponential and bimodal fits, overlays comparing conditions, and an adjustable histogram bin width. |
| 5 | [`05_spider_plots`](05_spider_plots) | Spider plots, with every trajectory starting at the origin: one per FOV, one combining all FOVs, a 2×2 comparison of conditions, and the top-N tracks by duration or path length. |

Each folder has its own README with inputs, settings and outputs.

## Installation

```bash
git clone https://github.com/diz-lab/Caging-of-membrane-to-cortex-attachment-proteins-can-trigger-cellular-symmetry-breaking-.git
cd Caging-of-membrane-to-cortex-attachment-proteins-can-trigger-cellular-symmetry-breaking-/smt-first-passage-analysis
pip install -r requirements.txt
```

Tested with Python 3.10, numpy 2.2, pandas 2.3, scipy 1.15, scikit-learn 1.7, matplotlib 3.10 and seaborn 0.13.

## Input data

The scripts read the CSV files exported from TrackMate:

- `FOV<n>-spots.csv`: one row per localization, with `TRACK_ID`, `POSITION_X`, `POSITION_Y`, `POSITION_T` and
  `FRAME`
- `FOV<n>-tracks.csv`: one row per track, with `TRACK_ID`, `TRACK_DURATION`, `TRACK_MEAN_SPEED` and
  `TRACK_DISPLACEMENT`

Every step expects one folder per condition:

```
data/
├── MoeS/        FOV1-spots.csv, FOV1-tracks.csv, …
├── MoeL/
├── MoeSActin/
└── MoeLActin/
```

Raw data is not included in this repository.

### Units and frame intervals

- **Step 1** expects the TrackMate export in **pixels and frames**. It converts them itself, with 0.097 µm/pixel
  and the frame interval you enter.
- **Step 2 onwards** expect positions in **µm** and times in **seconds**.
- **Frame intervals:** MoeL and MoeLActin were imaged at **20 ms** and MoeS and MoeSActin at **30 ms**. Run step 1
  separately for each frame interval. `sort_mobile_tracks.py` already uses 0.02 s and 0.03 s per condition.

## Usage

Every script is interactive: run it from a terminal and it asks for its folders and parameters. Press Enter to
accept a default value.

1. **Classify.** Run `01_trajectory_classification/trajectory_classification.py` on each condition folder of
   pixel/frame exports. The settings used for the paper were 0.097 µm/pixel, a frame interval of 0.02 s or 0.03 s,
   5 minimum localizations, and a minimum duration of 0.18 s. The minimum duration is inferred from the saved
   results; see the [step 1 README](01_trajectory_classification).
2. **Filter.** Put each `FOV<n>_comprehensive_results.csv` from step 1 in its condition folder, next to the
   µm/second spots and tracks files. Then run `02_median_displacement/sort_mobile_tracks.py` on the base folder and
   enter the duration threshold. The manuscript FPT/FET figures used 1.0 s; a 0.6 s set was also made. It writes
   `FOV<n>-0.2um-dur1s-mobile spots.csv` and `… tracks.csv` into each condition folder.
3. **Collect the filtered tracks.** Copy them to a new base folder with the same condition subfolders, and rename
   them to `FOV<n>-spots.csv` and `FOV<n>-tracks.csv`. Steps 3–5 only read files with these names.
4. **Analyse.** Run steps 3–5 on that new base folder. For the manuscript FPT/FET figures, set the bin width in
   step 4 as described in the [step 4 README](04_bimodal_fpt_fet), and enter the FET thresholds `0.5, 1.0`.

```bash
python 03_first_passage_time/first_passage_time.py
python 04_bimodal_fpt_fet/bimodal_fpt_fet_analysis.py
python 05_spider_plots/spider_plots.py
```

`02_median_displacement/median_displacement.py` is optional. It plots one condition folder at a time.

## Source files

| Repository file | Original file |
|---|---|
| `01_trajectory_classification/trajectory_classification.py` | `260120-UTA-correctedcode.py` |
| `02_median_displacement/sort_mobile_tracks.py` | `Sortingmobiletracks.docx` |
| `02_median_displacement/median_displacement.py` | `Median-displacement-gemini-code.docx` |
| `03_first_passage_time/first_passage_time.py` | `first-passage-timeanalysis-fixed.py` |
| `04_bimodal_fpt_fet/bimodal_fpt_fet_analysis.py` | `manuscript-complete_smt_analysis_v3.py` |
| `05_spider_plots/spider_plots.py` | `spider_plot_generator_albatalks.py` |

The code from the Word documents was copied as plain text without changes.

## Changes from the original code

Only one script was changed: step 4. Its histogram bin width and bin-count limits are now settings,
`BIN_WIDTH_FACTOR`, `MIN_BINS` and `MAX_BINS`. With their defaults, the script gives exactly the original results.
The manuscript figures were made with an edited copy that used other values, and the settings make that version
reproducible. All other scripts are unchanged.

## Checked against the saved results

The pipeline was rerun on part of the original data:

- **Step 1** reproduces the saved per-track values (MBCR, MBCC–CoM, MSD, diffusion coefficient, duration, speed,
  tortuosity) for all 814 tracks of MoeL FOV28. For MBCC–CoM, 4 tracks where two candidate circle centres tie
  depend on the order of the rows. The classification itself could not be compared, because only the Mobile tracks
  were available and the thresholds are fitted to all tracks in a FOV.
- **The filter** reproduces the saved track counts for every FOV of MoeL and MoeS. With a 1.0 s threshold, it selects
  exactly the tracks of the 1 s dataset used for the manuscript figures.
- **Step 4** reproduces all 20 statistics tables saved with the manuscript FPT/FET figures (`bins_half_width` and
  `bins_2x_width`), with the settings in its README.
- **Steps 3–5** reproduce the saved track counts for all four conditions.

## Known issues

- **Unit conversion between steps 1 and 2 is not part of this repository.** The files used for the paper are in µm
  and seconds, but `FRAME` holds times and `TRACK_X_LOCATION`/`TRACK_START` are still in pixels/frames, so they
  were not converted by `convert_units_fov.py`.
- **The filter's output file names** don't match the `FOV*-spots.csv` pattern that steps 3–5 look for, so they
  have to be renamed (usage step 3).
- **Floating-point edge case in `median_displacement.py`.** Tracks lasting exactly the minimum duration (e.g.
  0.6 s) can compute as 0.5999… s and be dropped. This does not affect the filter or steps 3–5.

## License

MIT, see [LICENSE](LICENSE).

## Author

Srishti Dar

Code development was supported by Claude (Anthropic).

