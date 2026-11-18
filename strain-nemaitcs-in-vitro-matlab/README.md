# MATLAB analysis code: actin networks on supported lipid bilayers (collaboration with Alba)

Generic, batch-ready versions of the MATLAB analyses used for this collaboration:
nematic order (alignment), PIV strain, temporal autocorrelation, bundling parameter and g(r).
Every analysis runs on any TIFF image or stack: set a file, a folder or a list at the top of the
`run_*.m` script, and each stack writes its own outputs next to it plus one summary CSV per folder.

## Quick start

```matlab
run('/path/to/Alba_MATLAB_analysis_code/setup_paths.m')   % once per session
edit run_PIV_strain          % set INPUT and parameters, then Run
```
Set `INPUT = ''` to get a folder picker. Pre-processing as in the Methods: ImageJ background
subtraction (rolling ball 50 px). For bundling, crop to 500 x 500 px. For g(r), binarise first.

## Folder layout

| Folder | Entry point | Engine (original code) | Per-stack outputs |
|---|---|---|---|
| `00_shared_utils` | - | `loadtiff.m`, `fminsearchbnd.m` (File Exchange), `llmFig.m` + new helpers `list_tiffs`, `read_tiff_stack`, `write_csv`, `set_defaults`, `save_figure` | - |
| `01_PIV_strain` | `run_PIV_strain.m` -> `piv_strain_analysis.m` | `src/` = lab cartPIV / MatPIV pipeline (63 files, unchanged) | `<name>_PIV/`: `ExpParams.mat`, `WS/cartPIV_WS.mat`, `PIV_strain.mat`, `strain_vs_time.csv/.pdf`, `strain_heatmap.pdf` |
| `02_alignment_nematic_order` | `run_alignment_nematic.m` -> `nematic_order_stack.m` | `FFTAlignment.m` (Oakes / Gardel lab, updated Linsmeier), `alignment.m` (Banerjee), unchanged | `<name>_alignment/`: `order_<t>.mat`, `align_<t>.mat`, `alignmentParams.mat`, `nematic_order.csv`, `nematicAlignment.pdf` |
| `03_temporal_autocorrelation` | `run_temporal_autocorr.m` -> `temporal_autocorr.m`; `plot_autocorr_overlay.m` | `ZGS_temp_autocorr.m` | `<name>_autocorr.mat/.csv/.pdf` |
| `04_bundling_parameter` | `run_bundling_parameter.m` -> `bundling_parameter.m` | two `threshold_code_ZGS.m` variants | `<name>_bundling.csv` |
| `05_pair_correlation_gofr` | `run_pair_correlation_gofr.m` -> `pair_correlation_gofr.m`; `plot_gofr_overlay.m` | `GrForXlinksAfines.m` (V. Yadav) | `<name>_gofr.mat/.csv/.pdf` |

Each analysis folder also has `original_scripts/`, holding the untouched historical scripts for
reference. `setup_paths.m` keeps them off the path.

## Methods and parameters

**01 PIV strain.** MatPIV single pass between consecutive frames (window 32 px, overlap 0.5),
then SNR (1.01), global (3 SD) and local median (3x3, 5 SD) filters and NaN interpolation.
Displacements are converted to um. The strain map is the divergence of the cumulative
displacement field (`PIVstrain.m`). The reported strain is
eps(t) = -<div x>, sign-flipped so that contraction is positive; the max strain is max eps(t).
The strain rate is the slope of eps(t) over a fit window chosen by `P.rateMethod`:
- `'pct25_75'` (default): from where eps first reaches 25 % to where it first reaches 75 % of
  its maximum (lab script `runForStrainAnalysis.m`).
- `'maxslope'`: the steepest window of `P.fitWindow_s` seconds.
- `'manual'`: `P.fitRange_s`.

Required parameters: `P.dscl` (um/px) and `P.tscl` (s/frame); the Alba data used 0.1076 and 10.5.
This replaces the copy-`DefineParams.m`-next-to-every-stack workflow. The same `ExpParams.mat`
is still written, so `runCartPIV.m` and the other `src/` tools keep working on the output folder.

**02 Nematic order.** Each window is Gaussian-filtered and transformed with a 2D FFT. The local
orientation is orthogonal to the axis of least second moment of the masked |FFT|. The local order
is n = 2(<cos^2(theta_j - theta_i)> - 1/2) over the (2 st + 1)^2 neighbouring windows (st = 1 gives
a 3x3 kernel). The per-frame value is `abs(nanmean(nanmean(order)))`, as in `geT_Nematic_ZGS.m`;
the CSV also gives the plain mean over all valid windows. Defaults are the values stored in the
revision `alignmentParams.mat`: winsize 19, overlap 0.5, st 1, mask_method 2. `checkpoint`
defaults to 0; the revision used 7220. It is an absolute summed intensity per window, so rescale
it for other bit depths or window sizes (7220 / 19^2 = 20 counts per pixel).
Requires the Image Processing Toolbox (`fspecial`, `improfile`, `bwlabel`, `regionprops`) and
`nanmean`. Outputs now go to `<name>_alignment/` instead of `alignment/`, so several stacks in one
folder no longer overwrite each other.

**03 Temporal autocorrelation.** Each frame is z-scored. For every lag of dF frames, C(dF) is the
mean over pairs of <z_k z_(k+dF)>. C is normalised by its first point and fitted over the first
`fitFraction` (default 0.3) of the lags to A e^(t/a2) + (1-A) e^(t/a3) with a2, a3 < 0
(`fminsearchbnd`). The script reports tau = -a2 and -a3, and tau_long = the larger of the two.
R^2 is computed as in the original. The `.mat` file keeps the variable names of
`actinCorrelations2Exp.mat` (`param`, `allCors`, `stats`).

**04 Bundling parameter.** lambda = nu / mu of the thresholded image, computed for every frame.
`P.method = 'normalize_total'` (default; Methods text): zeros become NaN, the image is divided by
its total intensity, and negatives become NaN. `'subtract_mean'`: subtract 0.9 x mean, then
negatives become NaN.

**05 g(r).** The foreground pixels of a binary image are treated as particles. For each interior
particle, g_i(r) = n(r, r+dr) / ring area / local density, where the local density is counted in a
(2 rlim)^2 box. g(r) is the mean of g_i(r) over particles, with the last 4 bins dropped.
g/g_max over bins 2-14 is fitted to A e^(-r/xi) + c. Defaults: rlim 35 px, dr 1 px,
`pixPerUm` 9.2937.

## Places where the original code and the Methods text differ

The defaults reproduce the original code, so earlier numbers can be regenerated. Options select
the alternative.

1. **Bundling, `statMode`.** The original computes `nu = nanstd(nanstd(I))`, which is the std across
   columns of the per-column std, and `mu = nanmean(nanmean(I))`, the mean of the column means.
   The Methods text describes the std and mean of the whole image. Set `P.statMode = 'global'` for
   the latter. Normalising by total intensity does not change lambda.
2. **g(r), `edgeMode`.** The original "interior particle" test, `-dim/2+rlim < x < dim/2-rlim`, was
   written for centred simulation coordinates. On image pixel coordinates (1..dim) it keeps the
   top-left quarter of the image, including pixels at the edge. `P.edgeMode = 'interior'` uses
   `rlim < x <= dim - rlim` on both axes. On a 120x140 test image the original rule kept 790
   particles and the interior rule kept 3,329.
3. **Autocorrelation.** The Methods equation has an offset a(4); the code has none. Set
   `P.fitOffset = true` to add it. The Methods mention fitting about 40 % of the curve; the code
   uses 30 % (`P.fitFraction`). C is normalised by the lag-1 point, not lag 0.
4. **PIV.** `runForStrainAnalysis.m` reported the max strain over frames 1-25 only. Here it is the
   maximum over the whole trace.

## Implementation notes

- TIFFs are read with `imfinfo`/`imread`, so any bit depth works. RGB images use their first channel.
- The g(r) particle search reads each particle's box directly from the mask instead of searching
  the full particle list. A Python check of both algorithms gave identical g(r) (max difference 0)
  in both edge modes on a random binary image.
- New code uses no toolbox functions, except the alignment engine (Image Processing Toolbox)
  and `nanmean` in the original engines.
- Static check: all 95 .m files parse with MISS_HIT `mh_lint`. The code has not yet been run in
  MATLAB; please test it on one known stack per analysis first.

## Provenance of the original files

| File | Found in (Dropbox unless noted) |
|---|---|
| PIV `src/*.m`, MatPIV `mask.m` | Desktop/PostDoc/Heisenberg_group/2026/PIV_code_ZGS/PIV_code/ |
| `DefineParams_Alba.m` | Collaboration_Alba/zach_analysis/contraction_analysis/DefineParams.m |
| `runForStrainAnalysis.m` | ForminProject/PIV Analysis/ |
| `alignment.m`, `FFTAlignment.m`, `geT_Nematic_ZGS.m`, `barPlotsAlignment_ZGS.m` | Manuscript_MyosinAssembly/Fig1_ActinGrowsFromSurface/prepfiles/Alignment_Zach/ |
| `ZGS_temp_autocorr.m` | Collaboration_Alba/zach_analysis/autocorrelation_analysis/ |
| `multi_plt_ZGS.m` | Collaboration_Alba/Thermal fluctuation/0Chol/E-21/Mov-1/Zach_analysis/ |
| `threshold_code_ZGS_subtractMean.m` | Collaboration_Alba/zach_analysis/image_bundle_threshold_ZGS/ |
| `threshold_code_ZGS_normalizeTotal.m` | Collaboration_Alba/zach_analysis/spatial_linescan_analysis/image_threshold_method/ |
| `GrForXlinksAfines.m`, `gofr_plt.m` | Collaboration_Alba/zach_analysis/g(r)_analysis/(code/) |
| `loadtiff.m` | Collaboration_Alba/zach_analysis/ (identical copy in g(r)_analysis/code) |
| `fminsearchbnd.m` | Manuscript_MyosinAssembly/code/ (identical in Manuscript_Rupture/code) |
| `llmFig.m` | Manuscript_Rupture/code/ |
