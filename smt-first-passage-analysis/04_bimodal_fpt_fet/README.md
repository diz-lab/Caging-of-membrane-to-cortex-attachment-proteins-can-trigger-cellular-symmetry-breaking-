# 4 · FPT / FET fitting (manuscript figures)

`bimodal_fpt_fet_analysis.py` extends step 3 with more radii and thresholds and an automatic choice of fit model.
It makes the first passage time and free exit time figures of the manuscript. It is
`manuscript-complete_smt_analysis_v3.py`, plus a bin-width setting (see below).

## Input

The same as step 3: a base folder containing `MoeS/`, `MoeL/`, `MoeSActin/` and `MoeLActin/`, each with the filtered
`FOV*-spots.csv` files. The manuscript figures used the tracks lasting at least **1 s** (the filter in step 2 with a
1.0 s threshold).

## Settings

**Asked when the script starts:**

- Base directory and output directory
- FET thresholds, comma-separated. The default is `0.3, 0.5, 1.0` µm, and 0.5 µm is always included. The manuscript
  figures used `0.5, 1.0`.

**Set at the top of the file:**

| Setting | Default | Manuscript figures |
|---|---|---|
| `BIN_WIDTH_FACTOR` | 1.0 | 0.5 for `bins_half_width`, 2 for `bins_2x_width` |
| `MIN_BINS`, `MAX_BINS` | 10, 100 | 5, 200 |
| `FPT_XLIM`, `FPT_YLIM`, `FET_XLIM`, `FET_YLIM` | axis limits of the 1 µm FPT and 0.5 µm FET figures | unchanged |

The FPT radii are fixed at 0.2, 0.6, 1, 2 and 3 µm.

## Method

The FPT and FET definitions are the same as in step 3. Each histogram is binned with the Freedman–Diaconis rule:
bin width 2·IQR/n^⅓ × `BIN_WIDTH_FACTOR`, giving between `MIN_BINS` and `MAX_BINS` bins. It is then fitted with:

- **Gaussian**
- **Exponential decay**
- **Bimodal (double Gaussian):** tried only when the smoothed histogram has two peaks, with the smaller at least 30%
  of the larger, separated by a dip of at least 30%. The reported mean and SD are the amplitude-weighted averages of
  the two peaks.

The fit with the highest R² is kept.

## Output

- `FPT_distribution_radius_<r>um.svg` and `FPT_statistics_radius_<r>um.csv`, for all radii
- `FET_distribution_threshold_<t>um.svg` and `FET_statistics_threshold_<t>um.csv`, for all thresholds
- `FPT_1.0um_individual.svg`, `_overlay.svg`, `_overlay_fits.svg` and `_statistics.csv`, on fixed axes
- `FET_0.5um_individual.svg`, `_overlay.svg`, `_overlay_fits.svg` and `_statistics.csv`, on fixed axes
- `end_to_end_distributions.svg`, `path_length_distributions.svg` and `lifetime_distributions.svg`
- `summary_statistics.csv`

The overlays plot MoeL over MoeLActin and MoeS over MoeSActin.

## Reproducing the manuscript figures

Run on the 1 s dataset with FET thresholds `0.5, 1.0`:

- With `BIN_WIDTH_FACTOR = 0.5`, `MIN_BINS = 5` and `MAX_BINS = 200`, the script reproduces all 10 statistics tables
  in `Manuscript-graphs-FPT-FET/bins_half_width`.
- With `BIN_WIDTH_FACTOR = 2` and the same limits, it reproduces all 10 tables in `bins_2x_width`.

Fit types, labels and track counts are identical. The fitted values differ by less than 10⁻⁷, except one SD at 0.2 µm
in `bins_2x_width`, which differs by 0.15% because the curve-fit converges slightly differently.

## Changes from the original code

The bin width and the bin-count limits were hard-coded (×1, 10–100 bins). They are now the three settings above, and
the defaults give exactly the original results. The manuscript figures were made with an edited copy of the script
that used other values; the settings make that version reproducible.
