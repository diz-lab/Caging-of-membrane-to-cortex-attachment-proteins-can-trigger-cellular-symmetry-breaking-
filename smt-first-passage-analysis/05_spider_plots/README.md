# 5 · Spider plots

`spider_plots.py` draws trajectories with every track moved so that it starts at the origin, for each condition.

## Input

The same as steps 3–4: a base folder containing `MoeS/`, `MoeL/`, `MoeSActin/` and `MoeLActin/`, each with
`FOV*-spots.csv` files.

## Settings (asked when the script starts)

| Prompt | Default |
|---|---|
| Base directory and output directory | none |
| Axis limit (µm) | 3.5 |
| Number of top tracks | 10 |
| Plot choice | 1–6, or comma-separated |

## Plot options

1. Spider plot for each FOV in each condition
2. Combined spider plot for each condition, with all FOVs merged
3. 2×2 comparison, with one panel per condition
4. Top-N longest tracks by duration (number of frames), one plot per track
5. Top-N longest tracks by path length, one plot per track
6. All of the above

## Output

- `<condition>_spider_plots/<FOV>_spider.png`
- `<condition>_COMBINED_spider.png`
- `ALL_CONDITIONS_comparison_spider.png`
- `<condition>_top<N>_duration/rank<k>_duration_<frames>frames.png`
- `<condition>_top<N>_pathlen/rank<k>_pathlen_<length>um.png`
