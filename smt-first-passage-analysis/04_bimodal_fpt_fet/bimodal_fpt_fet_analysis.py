"""
COMPLETE SMT ANALYSIS  v3
=========================

Fitting:      Gaussian vs Exponential from Revisedbins-gausvsexpofit.py
              + stringent bimodal detection (only attempted when dip test p<0.05
                AND a clear secondary local maximum exists above 30% of main peak)
Colours:      MoeL  #FF6347  MoeS  #FA8072
              MoeLActin #8B0000  MoeSActin #8B0000
Axes:         FPT@1µm  x=0–1.8s  y=0–3 PD  |  FET@0.5µm  x=0–8s  y=0–1.8 PD
              All FPT plots per radius: shared X and Y across conditions
              All FET plots per threshold: shared X and Y across conditions
Overlays:     MoeL over MoeLActin, MoeS over MoeSActin — histograms only, no fits
Fit stats:    printed as figure subtitle BELOW each subplot (not on the plot face)
Frame:        full four-sided frame on all plots, no grid inside
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from glob import glob
from scipy.optimize import curve_fit
import warnings
warnings.filterwarnings('ignore')

# ─────────────────────────────────────────────────────────────────────────────
# Colours
# ─────────────────────────────────────────────────────────────────────────────
CONDITION_COLORS = {
    'MoeL':      '#FF6347',   # tomato
    'MoeS':      '#FA8072',   # salmon
    'MoeLActin': '#8B0000',   # dark red
    'MoeSActin': '#8B0000',   # dark red  (same as MoeLActin per request)
}

OVERLAY_PAIRS = [
    ('MoeL', 'MoeLActin'),
    ('MoeS', 'MoeSActin'),
]

# Fixed axis limits
FET_SPECIAL_THRESHOLD = 0.5
FET_XLIM = (0, 8)     # seconds
FET_YLIM = (0, 1.8)

FPT_SPECIAL_RADIUS = 1.0
FPT_XLIM = (0, 1.8)   # seconds
FPT_YLIM = (0, 3.0)

# Histogram binning (Freedman-Diaconis width × BIN_WIDTH_FACTOR, limited to MIN_BINS–MAX_BINS)
# Defaults reproduce the original script. Manuscript figures: factor 0.5 (bins_half_width)
# or 2 (bins_2x_width), with MIN_BINS = 5 and MAX_BINS = 200.
BIN_WIDTH_FACTOR = 1.0
MIN_BINS = 10
MAX_BINS = 100


# ─────────────────────────────────────────────────────────────────────────────
# Fit functions  (from Revisedbins-gausvsexpofit.py, unchanged)
# ─────────────────────────────────────────────────────────────────────────────

def gaussian(x, amplitude, mean, stddev):
    return amplitude * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))


def double_gaussian(x, a1, m1, s1, a2, m2, s2):
    return gaussian(x, a1, m1, s1) + gaussian(x, a2, m2, s2)


def exponential_decay(x, amplitude, decay_rate, offset=0):
    return amplitude * np.exp(-decay_rate * x) + offset


def calculate_r_squared(y_actual, y_predicted):
    ss_res = np.sum((y_actual - y_predicted) ** 2)
    ss_tot = np.sum((y_actual - np.mean(y_actual)) ** 2)
    return 1 - (ss_res / ss_tot) if ss_tot != 0 else 0


def calculate_optimal_bins(data):
    """Freedman-Diaconis rule (from Revisedbins-gausvsexpofit.py)"""
    n = len(data)
    if n < 2:
        return MIN_BINS
    iqr = np.percentile(data, 75) - np.percentile(data, 25)
    if iqr == 0:
        return MIN_BINS
    bin_width = 2 * iqr / (n ** (1/3)) * BIN_WIDTH_FACTOR
    if bin_width <= 0:
        return MIN_BINS
    data_range = data.max() - data.min()
    n_bins = int(np.ceil(data_range / bin_width))
    return max(MIN_BINS, min(n_bins, MAX_BINS))


# ─────────────────────────────────────────────────────────────────────────────
# Stringent bimodal detection
# ─────────────────────────────────────────────────────────────────────────────

def _bimodal_guess(data, bin_centers, counts):
    """
    Return (guess, True) only if there is strong evidence of bimodality:
      - at least 2 local maxima in the smoothed histogram separated by a dip
        where the valley is < 70% of the SMALLER peak (i.e. dip is ≥ 30% deep)
      - the secondary peak is ≥ 30% of the primary peak height
    Returns (None, False) otherwise.
    """
    from scipy.signal import argrelextrema
    smooth = np.convolve(counts, np.ones(5) / 5, mode='same')
    peaks = argrelextrema(smooth, np.greater, order=3)[0]

    if len(peaks) < 2:
        return None, False

    # Keep only peaks >= 30 % of the tallest peak
    main_h = smooth[peaks].max()
    peaks = peaks[smooth[peaks] >= 0.30 * main_h]
    if len(peaks) < 2:
        return None, False

    # Sort by height, take the top two
    top2 = peaks[np.argsort(smooth[peaks])[-2:]]
    i1, i2 = sorted(top2)

    # Valley between the two peaks
    valley_h = smooth[i1:i2 + 1].min()
    smaller_peak_h = min(smooth[i1], smooth[i2])
    if valley_h > 0.70 * smaller_peak_h:   # dip not deep enough
        return None, False

    a1, m1 = float(smooth[i1]), float(bin_centers[i1])
    a2, m2 = float(smooth[i2]), float(bin_centers[i2])
    s_est = float(np.std(data)) * 0.4
    guess = [a1, m1, s_est, a2, m2, s_est]
    return guess, True


def fit_distribution(data, bin_centers, counts):
    """
    Try Gaussian and Exponential (from Revisedbins-gausvsexpofit.py).
    Additionally try bimodal only when stringent criteria are met.
    Return best result by R².
    """
    results = {}
    x_fine = np.linspace(data.min(), data.max(), 300)

    # ── Gaussian ──────────────────────────────────────────────────────────
    try:
        p0 = [counts.max(), float(np.mean(data)), float(np.std(data))]
        popt, _ = curve_fit(
            gaussian, bin_centers, counts, p0=p0, maxfev=10000,
            bounds=([0, data.min(), 0], [np.inf, data.max(), data.max()])
        )
        r2 = calculate_r_squared(counts, gaussian(bin_centers, *popt))
        results['gaussian'] = {
            'params': popt, 'r_squared': r2,
            'x_fit': x_fine, 'y_fit': gaussian(x_fine, *popt),
            'mean': popt[1], 'std': abs(popt[2]),
            'fit_type': 'Gaussian',
            'fit_color': 'black',
        }
    except Exception:
        pass

    # ── Exponential ───────────────────────────────────────────────────────
    try:
        decay_guess = 1.0 / np.mean(data) if np.mean(data) > 0 else 1.0
        p0 = [counts.max(), decay_guess, 0]
        popt, _ = curve_fit(
            exponential_decay, bin_centers, counts, p0=p0, maxfev=10000,
            bounds=([0, 0, 0], [np.inf, np.inf, counts.max()])
        )
        r2 = calculate_r_squared(counts, exponential_decay(bin_centers, *popt))
        exp_mean = 1.0 / popt[1] if popt[1] > 0 else np.nan
        results['exponential'] = {
            'params': popt, 'r_squared': r2,
            'x_fit': x_fine, 'y_fit': exponential_decay(x_fine, *popt),
            'mean': exp_mean, 'std': exp_mean,
            'fit_type': 'Exponential',
            'fit_color': 'black',
        }
    except Exception:
        pass

    # ── Bimodal  (stringent gate) ──────────────────────────────────────────
    guess, is_bimodal = _bimodal_guess(data, bin_centers, counts)
    if is_bimodal:
        try:
            lo = [0, data.min(), 0, 0, data.min(), 0]
            hi = [np.inf, data.max(), data.max(),
                  np.inf, data.max(), data.max()]
            popt, _ = curve_fit(double_gaussian, bin_centers, counts,
                                p0=guess, maxfev=15000, bounds=(lo, hi))
            r2 = calculate_r_squared(counts, double_gaussian(bin_centers, *popt))
            w_mean = ((popt[0] * popt[1] + popt[3] * popt[4])
                      / (popt[0] + popt[3] + 1e-12))
            w_std  = ((popt[0] * abs(popt[2]) + popt[3] * abs(popt[5]))
                      / (popt[0] + popt[3] + 1e-12))
            results['bimodal'] = {
                'params': popt, 'r_squared': r2,
                'x_fit': x_fine, 'y_fit': double_gaussian(x_fine, *popt),
                'mean': w_mean, 'std': w_std,
                'fit_type': 'Bimodal',
                'fit_color': 'black',
                'peaks': [(popt[1], abs(popt[2])), (popt[4], abs(popt[5]))],
            }
        except Exception:
            pass

    if not results:
        return None

    best_key = max(results, key=lambda k: results[k]['r_squared'])
    return results[best_key]


# ─────────────────────────────────────────────────────────────────────────────
# Data pipeline  (unchanged from v2 / Revisedbins script)
# ─────────────────────────────────────────────────────────────────────────────

def calculate_first_passage_time(track_data, cols, target_radius):
    x0 = track_data[cols['x']].iloc[0]
    y0 = track_data[cols['y']].iloc[0]
    t0 = track_data[cols['t']].iloc[0]
    distances = np.sqrt((track_data[cols['x']] - x0)**2 +
                        (track_data[cols['y']] - y0)**2)
    exit_points = distances > target_radius
    if exit_points.any():
        return track_data.loc[exit_points.idxmax(), cols['t']] - t0
    return np.nan


def calculate_free_exit_time(track_data, cols, confinement_threshold=0.5):
    displacements = np.sqrt(track_data[cols['x']].diff()**2 +
                            track_data[cols['y']].diff()**2)
    window_size = 3
    if len(displacements) < window_size:
        return np.nan
    rolling_mean = displacements.rolling(window=window_size).mean()
    transition_idx = (rolling_mean > confinement_threshold).idxmax()
    if pd.isna(transition_idx) or transition_idx == displacements.index[0]:
        return track_data[cols['t']].iloc[-1] - track_data[cols['t']].iloc[0]
    exit_time = track_data.loc[transition_idx, cols['t']] - track_data[cols['t']].iloc[0]
    return exit_time if exit_time > 0 else np.nan


def process_condition(condition_path, condition_name, radii, fet_thresholds):
    print(f"\n{'='*60}\nProcessing condition: {condition_name}\n{'='*60}")
    spots_files = sorted(glob(os.path.join(condition_path, "FOV*-spots.csv")))
    if not spots_files:
        print(f"WARNING: No spots files found in {condition_path}")
        return None
    print(f"Found {len(spots_files)} FOV files")

    all_fpt_data = {r: [] for r in radii}
    all_fet_data  = {t: [] for t in fet_thresholds}
    all_lifetimes, all_end_to_end, all_path_lengths = [], [], []

    for spots_file in spots_files:
        fov_id = os.path.basename(spots_file).split('-')[0]
        print(f"  Processing {fov_id}...")
        try:
            df = pd.read_csv(spots_file)
            cols = {
                'x':  'POSITION_X' if 'POSITION_X' in df.columns else 'x',
                'y':  'POSITION_Y' if 'POSITION_Y' in df.columns else 'y',
                't':  'POSITION_T' if 'POSITION_T' in df.columns else 't',
                'id': 'TRACK_ID'   if 'TRACK_ID'   in df.columns else 'track_id',
            }
            df = df.sort_values([cols['id'], cols['t']])

            def reset_time(g):
                g = g.copy(); g[cols['t']] -= g[cols['t']].iloc[0]; return g

            df = df.groupby(cols['id'], group_keys=False).apply(reset_time)

            def monotonic(g):
                t = g[cols['t']].values
                return all(t[i] <= t[i+1] for i in range(len(t)-1))

            df = df.groupby(cols['id'], group_keys=False).filter(monotonic)

            for _, td in df.groupby(cols['id']):
                if len(td) < 3:
                    continue
                for r in radii:
                    v = calculate_first_passage_time(td, cols, r)
                    if not np.isnan(v) and v > 0:
                        all_fpt_data[r].append(v)
                for t in fet_thresholds:
                    v = calculate_free_exit_time(td, cols, t)
                    if not np.isnan(v) and v > 0:
                        all_fet_data[t].append(v)
                all_lifetimes.append(td[cols['t']].iloc[-1] - td[cols['t']].iloc[0])
                x0, y0 = td[cols['x']].iloc[0], td[cols['y']].iloc[0]
                xf, yf = td[cols['x']].iloc[-1], td[cols['y']].iloc[-1]
                all_end_to_end.append(np.sqrt((xf-x0)**2 + (yf-y0)**2))
                dx, dy = td[cols['x']].diff(), td[cols['y']].diff()
                all_path_lengths.append(float(np.sqrt(dx**2+dy**2).sum()))
        except Exception as e:
            print(f"    Error {fov_id}: {e}")

    print(f"  Completed: {len(all_lifetimes)} tracks analysed")
    return {
        'condition_name': condition_name,
        'color': CONDITION_COLORS[condition_name],
        'fpt_data': all_fpt_data,
        'fet_data': all_fet_data,
        'lifetimes': all_lifetimes,
        'end_to_end': all_end_to_end,
        'path_lengths': all_path_lengths,
        'n_tracks': len(all_lifetimes),
    }


# ─────────────────────────────────────────────────────────────────────────────
# Plot helpers
# ─────────────────────────────────────────────────────────────────────────────

def _get(all_results, name):
    for r in all_results:
        if r is not None and r['condition_name'] == name:
            return r
    return None


def _apply_style(ax, xlabel, ylabel, title, xlim=None, ylim=None):
    """Common axis styling: full frame, no grid inside."""
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=13, fontweight='bold')
    ax.grid(False)
    # Keep all four spines so the plot has a complete frame
    for spine in ax.spines.values():
        spine.set_visible(True)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)


def _stats_caption(fit_result, n):
    """Build a one-line caption string to place BELOW the subplot."""
    if fit_result is None:
        return f"No fit converged   n = {n}"
    ft = fit_result['fit_type']
    r2 = fit_result['r_squared']
    mu = fit_result['mean']
    sd = fit_result['std']
    if ft == 'Bimodal' and 'peaks' in fit_result:
        p = fit_result['peaks']
        return (f"Fit: {ft}   Peak1: µ={p[0][0]:.3f} σ={p[0][1]:.3f}   "
                f"Peak2: µ={p[1][0]:.3f} σ={p[1][1]:.3f}   R²={r2:.3f}   n={n}")
    return f"Fit: {ft}   µ = {mu:.3f}   σ = {sd:.3f}   R² = {r2:.3f}   n = {n}"


def _plot_hist_fit(ax, data, color, n_bins, label='Data', alpha=0.65):
    """
    Draw histogram + best fit on ax.
    Returns (fit_result, caption_string).
    """
    counts, bins, _ = ax.hist(data, bins=n_bins, alpha=alpha, color=color,
                               edgecolor='black', linewidth=0.5,
                               density=True, label=label)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    fit_result = fit_distribution(data, bin_centers, counts)

    if fit_result:
        ft = fit_result['fit_type']
        lbl = f"{ft} fit (R²={fit_result['r_squared']:.3f})"
        ax.plot(fit_result['x_fit'], fit_result['y_fit'],
                color=fit_result['fit_color'], linewidth=2.5,
                linestyle='-', label=lbl)

    caption = _stats_caption(fit_result, len(data))
    return fit_result, caption


# ─────────────────────────────────────────────────────────────────────────────
# FPT @ 1 µm — individual 2×2 + overlays
# ─────────────────────────────────────────────────────────────────────────────

def plot_fpt_special(all_results, output_dir):
    radius = FPT_SPECIAL_RADIUS
    print(f"\nGenerating FPT special plots (radius = {radius} µm)...")
    cond_order = ['MoeL', 'MoeS', 'MoeLActin', 'MoeSActin']

    # ── Individual 2×2 ────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(14, 13))
    fig.suptitle(f'First Passage Time – Radius = {radius} µm  (Individual)',
                 fontsize=15, fontweight='bold')
    axes = axes.flatten()
    stats_summary = []

    for idx, cname in enumerate(cond_order):
        ax = axes[idx]
        result = _get(all_results, cname)
        if result is None:
            ax.text(0.5, 0.5, f'{cname}\nNot found', ha='center', va='center')
            ax.set_title(cname); continue

        data = np.array(result['fpt_data'].get(radius, []))
        data = data[np.isfinite(data)]
        data_vis = data[(data >= FPT_XLIM[0]) & (data <= FPT_XLIM[1])]

        if len(data_vis) < 3:
            ax.text(0.5, 0.5, f'Insufficient data\nn={len(data_vis)}',
                    ha='center', va='center')
            ax.set_title(cname); continue

        n_bins = calculate_optimal_bins(data_vis)
        fit_result, caption = _plot_hist_fit(ax, data_vis, result['color'], n_bins)
        _apply_style(ax, 'First Passage Time (s)', 'Probability Density',
                     cname, xlim=FPT_XLIM, ylim=FPT_YLIM)
        ax.legend(loc='upper right', fontsize=8, frameon=False)
        # Caption below subplot
        ax.set_xlabel(f'First Passage Time (s)\n{caption}', fontsize=9)

        if fit_result:
            stats_summary.append({
                'Condition': cname, 'Radius_um': radius,
                'Fit_Type': fit_result['fit_type'],
                'Mean_s': fit_result['mean'], 'Std_s': fit_result['std'],
                'R_squared': fit_result['r_squared'], 'N': len(data_vis),
            })

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'FPT_{radius}um_individual.svg'),
                format='svg', bbox_inches='tight')
    plt.close()

    # ── Overlay 1×2 (no fits, finer bins) ────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(f'First Passage Time – Radius = {radius} µm  (Overlay)',
                 fontsize=15, fontweight='bold')

    for ax, (light_name, dark_name) in zip(axes, OVERLAY_PAIRS):
        for cname, alpha, zorder in [(dark_name, 0.60, 2), (light_name, 0.60, 3)]:
            result = _get(all_results, cname)
            if result is None:
                continue
            data = np.array(result['fpt_data'].get(radius, []))
            data = data[np.isfinite(data)]
            data_vis = data[(data >= FPT_XLIM[0]) & (data <= FPT_XLIM[1])]
            if len(data_vis) < 3:
                continue
            n_bins = int(calculate_optimal_bins(data_vis) * 2)  # finer bins for overlay
            ax.hist(data_vis, bins=n_bins, alpha=alpha,
                    color=result['color'], edgecolor='black',
                    linewidth=0.5, density=True, label=cname, zorder=zorder)
        _apply_style(ax, 'First Passage Time (s)', 'Probability Density',
                     f'{light_name} vs {dark_name}',
                     xlim=FPT_XLIM, ylim=FPT_YLIM)
        ax.legend(loc='upper right', fontsize=9, frameon=False)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'FPT_{radius}um_overlay.svg'),
                format='svg', bbox_inches='tight')
    plt.close()

    # ── Overlay 1×2 WITH FITS ─────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(f'First Passage Time – Radius = {radius} µm  (Overlay + Fits)',
                 fontsize=15, fontweight='bold')

    _line_styles = ['-', '--']   # dark condition solid, light condition dashed

    for ax, (light_name, dark_name) in zip(axes, OVERLAY_PAIRS):
        for pair_idx, (cname, alpha, zorder) in enumerate(
                [(dark_name, 0.45, 2), (light_name, 0.45, 3)]):
            result = _get(all_results, cname)
            if result is None:
                continue
            data = np.array(result['fpt_data'].get(radius, []))
            data = data[np.isfinite(data)]
            data_vis = data[(data >= FPT_XLIM[0]) & (data <= FPT_XLIM[1])]
            if len(data_vis) < 3:
                continue
            n_bins = int(calculate_optimal_bins(data_vis) * 2)
            counts, bins, _ = ax.hist(
                data_vis, bins=n_bins, alpha=alpha,
                color=result['color'], edgecolor='black',
                linewidth=0.5, density=True, label=cname, zorder=zorder)
            bin_centers = (bins[:-1] + bins[1:]) / 2
            fit_result = fit_distribution(data_vis, bin_centers, counts)
            if fit_result:
                ls = _line_styles[pair_idx]
                lbl = f"{cname} {fit_result['fit_type']} (R²={fit_result['r_squared']:.3f})"
                ax.plot(fit_result['x_fit'], fit_result['y_fit'],
                        color=result['color'], linewidth=2.5,
                        linestyle=ls, label=lbl, zorder=zorder + 1)
        _apply_style(ax, 'First Passage Time (s)', 'Probability Density',
                     f'{light_name} vs {dark_name}',
                     xlim=FPT_XLIM, ylim=FPT_YLIM)
        ax.legend(loc='upper right', fontsize=8, frameon=False)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'FPT_{radius}um_overlay_fits.svg'),
                format='svg', bbox_inches='tight')
    plt.close()

    if stats_summary:
        pd.DataFrame(stats_summary).to_csv(
            os.path.join(output_dir, f'FPT_{radius}um_statistics.csv'), index=False)
    print(f"  FPT {radius} µm plots saved.")


# ─────────────────────────────────────────────────────────────────────────────
# FET @ 0.5 µm — individual 2×2 + overlays
# ─────────────────────────────────────────────────────────────────────────────

def plot_fet_special(all_results, output_dir):
    thr = FET_SPECIAL_THRESHOLD
    print(f"\nGenerating FET special plots (threshold = {thr} µm)...")
    cond_order = ['MoeL', 'MoeS', 'MoeLActin', 'MoeSActin']

    # ── Individual 2×2 ────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(14, 13))
    fig.suptitle(f'Free Exit Time – Threshold = {thr} µm  (Individual)',
                 fontsize=15, fontweight='bold')
    axes = axes.flatten()
    stats_summary = []

    for idx, cname in enumerate(cond_order):
        ax = axes[idx]
        result = _get(all_results, cname)
        if result is None:
            ax.text(0.5, 0.5, f'{cname}\nNot found', ha='center', va='center')
            ax.set_title(cname); continue

        data = np.array(result['fet_data'].get(thr, []))
        data = data[np.isfinite(data)]
        data_vis = data[(data >= FET_XLIM[0]) & (data <= FET_XLIM[1])]

        if len(data_vis) < 3:
            ax.text(0.5, 0.5, f'Insufficient data\nn={len(data_vis)}',
                    ha='center', va='center')
            ax.set_title(cname); continue

        n_bins = calculate_optimal_bins(data_vis)
        fit_result, caption = _plot_hist_fit(ax, data_vis, result['color'], n_bins)
        _apply_style(ax, 'Free Exit Time (s)', 'Probability Density',
                     cname, xlim=FET_XLIM, ylim=FET_YLIM)
        ax.legend(loc='upper right', fontsize=8, frameon=False)
        ax.set_xlabel(f'Free Exit Time (s)\n{caption}', fontsize=9)

        if fit_result:
            stats_summary.append({
                'Condition': cname, 'Threshold_um': thr,
                'Fit_Type': fit_result['fit_type'],
                'Mean_s': fit_result['mean'], 'Std_s': fit_result['std'],
                'R_squared': fit_result['r_squared'], 'N': len(data_vis),
            })

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'FET_{thr}um_individual.svg'),
                format='svg', bbox_inches='tight')
    plt.close()

    # ── Overlay 1×2 (no fits, finer bins) ────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(f'Free Exit Time – Threshold = {thr} µm  (Overlay)',
                 fontsize=15, fontweight='bold')

    for ax, (light_name, dark_name) in zip(axes, OVERLAY_PAIRS):
        for cname, alpha, zorder in [(dark_name, 0.60, 2), (light_name, 0.60, 3)]:
            result = _get(all_results, cname)
            if result is None:
                continue
            data = np.array(result['fet_data'].get(thr, []))
            data = data[np.isfinite(data)]
            data_vis = data[(data >= FET_XLIM[0]) & (data <= FET_XLIM[1])]
            if len(data_vis) < 3:
                continue
            n_bins = int(calculate_optimal_bins(data_vis) * 2)  # finer bins for overlay
            ax.hist(data_vis, bins=n_bins, alpha=alpha,
                    color=result['color'], edgecolor='black',
                    linewidth=0.5, density=True, label=cname, zorder=zorder)
        _apply_style(ax, 'Free Exit Time (s)', 'Probability Density',
                     f'{light_name} vs {dark_name}',
                     xlim=FET_XLIM, ylim=FET_YLIM)
        ax.legend(loc='upper right', fontsize=9, frameon=False)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'FET_{thr}um_overlay.svg'),
                format='svg', bbox_inches='tight')
    plt.close()

    # ── Overlay 1×2 WITH FITS ─────────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(f'Free Exit Time – Threshold = {thr} µm  (Overlay + Fits)',
                 fontsize=15, fontweight='bold')

    _line_styles = ['-', '--']   # dark condition solid, light condition dashed

    for ax, (light_name, dark_name) in zip(axes, OVERLAY_PAIRS):
        for pair_idx, (cname, alpha, zorder) in enumerate(
                [(dark_name, 0.45, 2), (light_name, 0.45, 3)]):
            result = _get(all_results, cname)
            if result is None:
                continue
            data = np.array(result['fet_data'].get(thr, []))
            data = data[np.isfinite(data)]
            data_vis = data[(data >= FET_XLIM[0]) & (data <= FET_XLIM[1])]
            if len(data_vis) < 3:
                continue
            n_bins = int(calculate_optimal_bins(data_vis) * 2)
            counts, bins, _ = ax.hist(
                data_vis, bins=n_bins, alpha=alpha,
                color=result['color'], edgecolor='black',
                linewidth=0.5, density=True, label=cname, zorder=zorder)
            bin_centers = (bins[:-1] + bins[1:]) / 2
            fit_result = fit_distribution(data_vis, bin_centers, counts)
            if fit_result:
                ls = _line_styles[pair_idx]
                lbl = f"{cname} {fit_result['fit_type']} (R²={fit_result['r_squared']:.3f})"
                ax.plot(fit_result['x_fit'], fit_result['y_fit'],
                        color=result['color'], linewidth=2.5,
                        linestyle=ls, label=lbl, zorder=zorder + 1)
        _apply_style(ax, 'Free Exit Time (s)', 'Probability Density',
                     f'{light_name} vs {dark_name}',
                     xlim=FET_XLIM, ylim=FET_YLIM)
        ax.legend(loc='upper right', fontsize=8, frameon=False)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f'FET_{thr}um_overlay_fits.svg'),
                format='svg', bbox_inches='tight')
    plt.close()

    if stats_summary:
        pd.DataFrame(stats_summary).to_csv(
            os.path.join(output_dir, f'FET_{thr}um_statistics.csv'), index=False)
    print(f"  FET {thr} µm plots saved.")


# ─────────────────────────────────────────────────────────────────────────────
# All-radii FPT  &  all-thresholds FET  (standard 2×2, no grid)
# ─────────────────────────────────────────────────────────────────────────────

def plot_fpt_distributions(all_results, radii, output_dir):
    print("\nGenerating FPT distributions (all radii)...")
    for radius in radii:
        # ── Compute shared axis limits across all conditions ───────────────
        all_x, all_pd = [], []
        for result in all_results:
            if result is None:
                continue
            d = np.array(result['fpt_data'].get(radius, []))
            d = d[np.isfinite(d)]
            if len(d) < 3:
                continue
            all_x.append(np.percentile(d, 99))
            # estimate max probability density via histogram
            nb = calculate_optimal_bins(d)
            counts_tmp, _ = np.histogram(d, bins=nb, density=True)
            all_pd.append(counts_tmp.max())
        shared_xlim = (0, max(all_x) * 1.05) if all_x else None
        shared_ylim = (0, max(all_pd) * 1.15) if all_pd else None

        fig, axes = plt.subplots(2, 2, figsize=(14, 13))
        fig.suptitle(f'First Passage Time  (Radius = {radius} µm)',
                     fontsize=16, fontweight='bold')
        axes = axes.flatten()
        stats_summary = []

        for idx, result in enumerate(all_results):
            if result is None:
                continue
            ax = axes[idx]
            data = np.array(result['fpt_data'].get(radius, []))
            data = data[np.isfinite(data)]
            if len(data) < 3:
                ax.text(0.5, 0.5, f'Insufficient Data\nn={len(data)}',
                        ha='center', va='center', fontsize=14)
                ax.set_title(result['condition_name']); continue

            n_bins = calculate_optimal_bins(data)
            fit_result, caption = _plot_hist_fit(ax, data, result['color'], n_bins)
            _apply_style(ax, f'First Passage Time (s)\n{caption}',
                         'Probability Density', result['condition_name'],
                         xlim=shared_xlim, ylim=shared_ylim)
            ax.legend(loc='upper right', fontsize=8, frameon=False)

            if fit_result:
                stats_summary.append({
                    'Condition': result['condition_name'], 'Radius_um': radius,
                    'Fit_Type': fit_result['fit_type'],
                    'Mean_s': fit_result['mean'], 'StdDev_s': fit_result['std'],
                    'R_squared': fit_result['r_squared'],
                    'N_tracks': len(data), 'N_bins': n_bins,
                })

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'FPT_distribution_radius_{radius}um.svg'),
                    format='svg', bbox_inches='tight')
        plt.close()
        if stats_summary:
            pd.DataFrame(stats_summary).to_csv(
                os.path.join(output_dir, f'FPT_statistics_radius_{radius}um.csv'),
                index=False)


def plot_fet_distributions(all_results, fet_thresholds, output_dir):
    print("\nGenerating FET distributions (all thresholds)...")
    for threshold in fet_thresholds:
        # ── Compute shared axis limits across all conditions ───────────────
        all_x, all_pd = [], []
        for result in all_results:
            if result is None:
                continue
            d = np.array(result['fet_data'].get(threshold, []))
            d = d[np.isfinite(d)]
            if len(d) < 3:
                continue
            all_x.append(np.percentile(d, 99))
            nb = calculate_optimal_bins(d)
            counts_tmp, _ = np.histogram(d, bins=nb, density=True)
            all_pd.append(counts_tmp.max())
        shared_xlim = (0, max(all_x) * 1.05) if all_x else None
        shared_ylim = (0, max(all_pd) * 1.15) if all_pd else None

        fig, axes = plt.subplots(2, 2, figsize=(14, 13))
        fig.suptitle(f'Free Exit Time  (Threshold = {threshold} µm)',
                     fontsize=16, fontweight='bold')
        axes = axes.flatten()
        stats_summary = []

        for idx, result in enumerate(all_results):
            if result is None:
                continue
            ax = axes[idx]
            data = np.array(result['fet_data'].get(threshold, []))
            data = data[np.isfinite(data)]
            if len(data) < 3:
                ax.text(0.5, 0.5, f'Insufficient Data\nn={len(data)}',
                        ha='center', va='center', fontsize=14)
                ax.set_title(result['condition_name']); continue

            n_bins = calculate_optimal_bins(data)
            fit_result, caption = _plot_hist_fit(ax, data, result['color'], n_bins)
            _apply_style(ax, f'Free Exit Time (s)\n{caption}',
                         'Probability Density', result['condition_name'],
                         xlim=shared_xlim, ylim=shared_ylim)
            ax.legend(fontsize=8, frameon=False)

            if fit_result:
                stats_summary.append({
                    'Condition': result['condition_name'], 'Threshold_um': threshold,
                    'Fit_Type': fit_result['fit_type'],
                    'Mean_s': fit_result['mean'], 'StdDev_s': fit_result['std'],
                    'R_squared': fit_result['r_squared'], 'N_tracks': len(data),
                })

        plt.tight_layout()
        plt.savefig(os.path.join(output_dir,
                                  f'FET_distribution_threshold_{threshold}um.svg'),
                    format='svg', bbox_inches='tight')
        plt.close()
        if stats_summary:
            pd.DataFrame(stats_summary).to_csv(
                os.path.join(output_dir,
                             f'FET_statistics_threshold_{threshold}um.csv'),
                index=False)


# ─────────────────────────────────────────────────────────────────────────────
# Other distributions  (end-to-end, path length, lifetime)  — no grid
# ─────────────────────────────────────────────────────────────────────────────

def plot_end_to_end_distributions(all_results, output_dir):
    print("\nGenerating end-to-end distance distributions...")
    all_d = [v for r in all_results if r for v in r['end_to_end']]
    if not all_d:
        return
    bins = np.linspace(0, np.percentile(all_d, 99), 50)
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle('End-to-End Distance Distribution', fontsize=16, fontweight='bold')
    for idx, result in enumerate(all_results):
        if result is None:
            continue
        ax = axes.flatten()[idx]
        data = np.array(result['end_to_end'])
        ax.hist(data, bins=bins, alpha=0.7, color=result['color'],
                edgecolor='black', linewidth=0.5)
        mu, med = np.mean(data), np.median(data)
        ax.axvline(mu,  color='black', linestyle='--', linewidth=1.5,
                   label=f'Mean: {mu:.2f}µm')
        ax.axvline(med, color='grey',  linestyle=':',  linewidth=1.5,
                   label=f'Median: {med:.2f}µm')
        _apply_style(ax, 'End-to-End Distance (µm)', 'Count',
                     result['condition_name'])
        ax.legend(fontsize=8, frameon=False)
        caption = f"n = {len(data)}   µ = {mu:.2f} µm   σ = {np.std(data):.2f} µm"
        ax.set_xlabel(f'End-to-End Distance (µm)\n{caption}', fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'end_to_end_distributions.svg'),
                format='svg', bbox_inches='tight')
    plt.close()


def plot_path_length_distributions(all_results, output_dir):
    print("\nGenerating path length distributions...")
    all_l = [v for r in all_results if r for v in r['path_lengths']]
    if not all_l:
        return
    bins = np.linspace(0, np.percentile(all_l, 99), 50)
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle('Path Length Distribution', fontsize=16, fontweight='bold')
    for idx, result in enumerate(all_results):
        if result is None:
            continue
        ax = axes.flatten()[idx]
        data = np.array(result['path_lengths'])
        ax.hist(data, bins=bins, alpha=0.7, color=result['color'],
                edgecolor='black', linewidth=0.5)
        mu, med = np.mean(data), np.median(data)
        ax.axvline(mu,  color='black', linestyle='--', linewidth=1.5,
                   label=f'Mean: {mu:.2f}µm')
        ax.axvline(med, color='grey',  linestyle=':',  linewidth=1.5,
                   label=f'Median: {med:.2f}µm')
        _apply_style(ax, 'Path Length (µm)', 'Count', result['condition_name'])
        ax.legend(fontsize=8, frameon=False)
        caption = f"n = {len(data)}   µ = {mu:.2f} µm   σ = {np.std(data):.2f} µm"
        ax.set_xlabel(f'Path Length (µm)\n{caption}', fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'path_length_distributions.svg'),
                format='svg', bbox_inches='tight')
    plt.close()


def plot_lifetime_distributions(all_results, output_dir):
    print("\nGenerating lifetime distributions...")
    all_lt = [v for r in all_results if r for v in r['lifetimes']]
    if not all_lt:
        return
    bins = np.linspace(0, np.percentile(all_lt, 99), 50)
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle('Track Lifetime Distribution', fontsize=16, fontweight='bold')
    for idx, result in enumerate(all_results):
        if result is None:
            continue
        ax = axes.flatten()[idx]
        data = np.array(result['lifetimes'])
        ax.hist(data, bins=bins, alpha=0.7, color=result['color'],
                edgecolor='black', linewidth=0.5)
        mu, med = np.mean(data), np.median(data)
        ax.axvline(mu,  color='black', linestyle='--', linewidth=1.5,
                   label=f'Mean: {mu:.2f}s')
        ax.axvline(med, color='grey',  linestyle=':',  linewidth=1.5,
                   label=f'Median: {med:.2f}s')
        _apply_style(ax, 'Lifetime (s)', 'Count', result['condition_name'])
        ax.legend(fontsize=8, frameon=False)
        caption = f"n = {len(data)}   µ = {mu:.2f} s   σ = {np.std(data):.2f} s"
        ax.set_xlabel(f'Lifetime (s)\n{caption}', fontsize=9)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'lifetime_distributions.svg'),
                format='svg', bbox_inches='tight')
    plt.close()


# ─────────────────────────────────────────────────────────────────────────────
# Summary CSV
# ─────────────────────────────────────────────────────────────────────────────

def generate_summary_report(all_results, output_dir):
    print("\nGenerating summary report...")
    rows = []
    for result in all_results:
        if result is None:
            continue
        first_thr = list(result['fet_data'].keys())[0]
        fv = result['fet_data'][first_thr]
        rows.append({
            'Condition': result['condition_name'],
            'N_Tracks': result['n_tracks'],
            'Mean_Lifetime_s': np.mean(result['lifetimes']),
            'Std_Lifetime_s':  np.std(result['lifetimes']),
            'Mean_EndToEnd_um': np.mean(result['end_to_end']),
            'Std_EndToEnd_um':  np.std(result['end_to_end']),
            'Mean_PathLength_um': np.mean(result['path_lengths']),
            'Std_PathLength_um':  np.std(result['path_lengths']),
            'Mean_FreeExitTime_s': np.mean(fv) if fv else np.nan,
            'Std_FreeExitTime_s':  np.std(fv)  if fv else np.nan,
        })
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(output_dir, 'summary_statistics.csv'), index=False)
    print(f"Summary saved.\n{'='*60}")
    print(df.to_string(index=False))
    print('='*60)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    print("\n" + "="*60)
    print("COMPLETE SMT ANALYSIS v3")
    print("="*60)

    base_dir   = input("\nEnter base directory (condition folders): ").strip().replace('"','').replace("'",'')
    output_dir = input("Enter output directory: ").strip().replace('"','').replace("'",'')

    print("\nFET thresholds (default: 0.3, 0.5, 1.0 — press Enter):")
    t_in = input("Thresholds (comma-separated): ").strip()
    fet_thresholds = ([float(x.strip()) for x in t_in.split(',')]
                      if t_in else [0.3, 0.5, 1.0])
    if FET_SPECIAL_THRESHOLD not in fet_thresholds:
        fet_thresholds.append(FET_SPECIAL_THRESHOLD)
    print(f"FET thresholds: {fet_thresholds} µm")

    os.makedirs(output_dir, exist_ok=True)

    conditions = ['MoeS', 'MoeL', 'MoeSActin', 'MoeLActin']
    radii = [0.2, 0.6, 1, 2, 3]
    if FPT_SPECIAL_RADIUS not in radii:
        radii.append(FPT_SPECIAL_RADIUS)

    all_results = []
    for cname in conditions:
        cpath = os.path.join(base_dir, cname)
        if not os.path.exists(cpath):
            print(f"\nWARNING: {cpath} not found")
            all_results.append(None)
            continue
        all_results.append(process_condition(cpath, cname, radii, fet_thresholds))

    print("\n" + "="*60 + "\nGENERATING PLOTS\n" + "="*60)
    plot_fpt_distributions(all_results, radii, output_dir)
    plot_fet_distributions(all_results, fet_thresholds, output_dir)
    plot_fpt_special(all_results, output_dir)
    plot_fet_special(all_results, output_dir)
    plot_end_to_end_distributions(all_results, output_dir)
    plot_path_length_distributions(all_results, output_dir)
    plot_lifetime_distributions(all_results, output_dir)
    generate_summary_report(all_results, output_dir)

    print("\n" + "="*60)
    print("ANALYSIS COMPLETE!  Results saved to:", output_dir)
    print("="*60)


if __name__ == "__main__":
    main()
