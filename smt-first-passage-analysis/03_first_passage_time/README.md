# 3 · First passage time

`first_passage_time.py` computes first passage times and other track metrics for four conditions.

## Input

A base folder with one subfolder per condition: `MoeS`, `MoeL`, `MoeSActin` and `MoeLActin`. Each contains
`FOV*-spots.csv` files.

## Settings (asked when the script starts)

- Base directory with the condition folders
- Output directory

The radii (1, 2 and 3 µm), the condition names and colours, and the FET threshold (0.5 µm) are set in the code.

## Method

For every track with at least 3 points, with time reset to 0 at the first point and tracks with non-monotonic time
removed:

- **First passage time (FPT):** the time until the particle is first further than *r* from its starting position.
  Tracks that never leave the radius are not counted.
- **Free exit time (FET):** the time until the rolling mean of 3 consecutive step sizes first exceeds 0.5 µm. Tracks
  with no such transition get their full duration.
- **Track lifetime, end-to-end distance and path length.**

The FPT histograms are fitted with a single Gaussian.

## Output

- `FPT_distribution_radius_<r>um.png` and `FPT_statistics_radius_<r>um.csv`
- `lifetime_distributions.png`, `end_to_end_distributions.png`, `path_length_distributions.png` and
  `free_exit_time_distributions.png`
- `summary_statistics.csv`
