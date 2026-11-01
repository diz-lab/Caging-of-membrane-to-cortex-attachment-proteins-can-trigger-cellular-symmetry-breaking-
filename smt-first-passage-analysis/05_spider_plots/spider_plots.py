import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
from glob import glob
import warnings
warnings.filterwarnings('ignore')

# ──────────────────────────────────────────────────────────────────────────────
# Updated color / linewidth scheme
#   MoeLActin & MoeSActin → dark red,  lw = 2.4  (double of original 1.2)
#   MoeS      & MoeL      → red,       lw = 2.4
# ──────────────────────────────────────────────────────────────────────────────
CONDITIONS = {
    'MoeS':      {'color': 'red',     'lw': 2.4},
    'MoeL':      {'color': 'red',     'lw': 2.4},
    'MoeSActin': {'color': 'darkred', 'lw': 2.4},
    'MoeLActin': {'color': 'darkred', 'lw': 2.4},
}

# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

def detect_cols(df):
    return {
        'x':  'POSITION_X' if 'POSITION_X' in df.columns else 'x',
        'y':  'POSITION_Y' if 'POSITION_Y' in df.columns else 'y',
        't':  'POSITION_T' if 'POSITION_T' in df.columns else 't',
        'id': 'TRACK_ID'   if 'TRACK_ID'   in df.columns else 'track_id',
    }


def load_all_tracks(condition_path):
    """Return list of (x_norm, y_norm, n_frames, total_path_length) for every valid track."""
    spots_files = sorted(glob(os.path.join(condition_path, "FOV*-spots.csv")))
    records = []
    for f in spots_files:
        try:
            df = pd.read_csv(f)
            cols = detect_cols(df)
            df = df.sort_values([cols['id'], cols['t']])
            for _, td in df.groupby(cols['id']):
                if len(td) < 2:
                    continue
                x0, y0 = td[cols['x']].iloc[0], td[cols['y']].iloc[0]
                xn = td[cols['x']].values - x0
                yn = td[cols['y']].values - y0
                # cumulative path length (arc length)
                dx = np.diff(xn)
                dy = np.diff(yn)
                path_len = float(np.sum(np.sqrt(dx**2 + dy**2)))
                records.append({
                    'xn': xn, 'yn': yn,
                    'n_frames': len(td),
                    'path_len': path_len,
                })
        except Exception as e:
            print(f"    Warning – could not load {f}: {e}")
    return records, spots_files


def _add_reference_circles(ax, axis_limit, radii=(1, 2, 3)):
    """Add dotted reference circles and cross-hairs."""
    ax.axhline(0, color='black', linewidth=1.5, alpha=0.5, zorder=0)
    ax.axvline(0, color='black', linewidth=1.5, alpha=0.5, zorder=0)
    ax.grid(True, linestyle=':', alpha=0.4, color='gray', linewidth=0.8)
    for r in radii:
        if r <= axis_limit:
            circle = plt.Circle((0, 0), r, color='gray', fill=False,
                                 linestyle=':', linewidth=1.2, alpha=0.5, zorder=0)
            ax.add_patch(circle)
            ax.text(r, 0.08, f'{r}µm', fontsize=9, color='gray',
                    ha='center', va='bottom')


def _style_ax(ax, axis_limit):
    ax.set_xlim(-axis_limit, axis_limit)
    ax.set_ylim(-axis_limit, axis_limit)
    ax.set_aspect('equal')
    ax.set_xlabel('ΔX (µm)', fontsize=13, fontweight='bold')
    ax.set_ylabel('ΔY (µm)', fontsize=13, fontweight='bold')


def _single_track_spider(ax, xn, yn, color, lw, axis_limit, label=''):
    _add_reference_circles(ax, axis_limit)
    ax.plot(xn, yn, color=color, linewidth=lw, alpha=0.85, zorder=2)
    ax.plot(0, 0, 'o', color='green', markersize=6, alpha=0.7, zorder=4, label='Start')
    ax.plot(xn[-1], yn[-1], 'o', color=color, markersize=5, alpha=0.8, zorder=5, label='End')
    _style_ax(ax, axis_limit)
    if label:
        ax.set_title(label, fontsize=11, fontweight='bold', pad=8)


# ──────────────────────────────────────────────────────────────────────────────
# NEW: Plot top-10 longest tracks by DURATION (n_frames) – one plot per track
# ──────────────────────────────────────────────────────────────────────────────

def plot_top10_longest_duration(condition_path, condition_name, color, lw,
                                 output_dir, axis_limit=3.5, n_top=10):
    """
    For a given condition, pick the top-n_top tracks with the most frames
    and save each as an individual spider plot.
    """
    print(f"  [Duration top-{n_top}] Processing {condition_name}...")
    records, spots_files = load_all_tracks(condition_path)
    if not records:
        print(f"    No tracks found – skipping.")
        return

    # Sort descending by n_frames
    records.sort(key=lambda r: r['n_frames'], reverse=True)
    top_tracks = records[:n_top]

    sub_dir = os.path.join(output_dir, f"{condition_name}_top{n_top}_duration")
    os.makedirs(sub_dir, exist_ok=True)

    for rank, rec in enumerate(top_tracks, start=1):
        fig, ax = plt.subplots(figsize=(8, 8))
        title = (f"{condition_name} | Top-{rank} by Duration\n"
                 f"Frames = {rec['n_frames']}  |  Path = {rec['path_len']:.2f} µm")
        _single_track_spider(ax, rec['xn'], rec['yn'], color, lw, axis_limit, label=title)
        ax.legend(loc='upper right', fontsize=9, framealpha=0.9)
        plt.tight_layout()
        out = os.path.join(sub_dir, f"rank{rank:02d}_duration_{rec['n_frames']}frames.png")
        plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        print(f"    Rank {rank:2d}: {rec['n_frames']} frames  → {os.path.basename(out)}")

    print(f"    Saved to {sub_dir}")


# ──────────────────────────────────────────────────────────────────────────────
# NEW: Plot top-10 longest tracks by PATH LENGTH – one plot per track
# ──────────────────────────────────────────────────────────────────────────────

def plot_top10_longest_path(condition_path, condition_name, color, lw,
                             output_dir, axis_limit=3.5, n_top=10):
    """
    For a given condition, pick the top-n_top tracks with the greatest
    total arc-length (cumulative step distances) and save each as an
    individual spider plot.
    """
    print(f"  [Path-length top-{n_top}] Processing {condition_name}...")
    records, _ = load_all_tracks(condition_path)
    if not records:
        print(f"    No tracks found – skipping.")
        return

    # Sort descending by path_len
    records.sort(key=lambda r: r['path_len'], reverse=True)
    top_tracks = records[:n_top]

    sub_dir = os.path.join(output_dir, f"{condition_name}_top{n_top}_pathlen")
    os.makedirs(sub_dir, exist_ok=True)

    for rank, rec in enumerate(top_tracks, start=1):
        fig, ax = plt.subplots(figsize=(8, 8))
        title = (f"{condition_name} | Top-{rank} by Path Length\n"
                 f"Path = {rec['path_len']:.2f} µm  |  Frames = {rec['n_frames']}")
        _single_track_spider(ax, rec['xn'], rec['yn'], color, lw, axis_limit, label=title)
        ax.legend(loc='upper right', fontsize=9, framealpha=0.9)
        plt.tight_layout()
        out = os.path.join(sub_dir, f"rank{rank:02d}_pathlen_{rec['path_len']:.2f}um.png")
        plt.savefig(out, dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        print(f"    Rank {rank:2d}: {rec['path_len']:.3f} µm  → {os.path.basename(out)}")

    print(f"    Saved to {sub_dir}")


# ──────────────────────────────────────────────────────────────────────────────
# Original functions (colours/lw updated via parameters)
# ──────────────────────────────────────────────────────────────────────────────

def generate_spider_plots_for_condition(condition_path, condition_name, color, lw,
                                         output_dir, axis_limit=3.5):
    print(f"\n{'='*60}")
    print(f"Generating spider plots for: {condition_name}")
    print(f"{'='*60}")

    spots_files = sorted(glob(os.path.join(condition_path, "FOV*-spots.csv")))
    if not spots_files:
        print(f"WARNING: No spots files found in {condition_path}")
        return
    print(f"Found {len(spots_files)} FOV files")

    condition_output_dir = os.path.join(output_dir, f"{condition_name}_spider_plots")
    os.makedirs(condition_output_dir, exist_ok=True)

    for spots_file in spots_files:
        fov_id = os.path.basename(spots_file).split('-')[0]
        print(f"  Processing {fov_id}...")
        try:
            df_spots = pd.read_csv(spots_file)
            cols = detect_cols(df_spots)
            df_spots = df_spots.sort_values([cols['id'], cols['t']])
            create_spider_plot(df_spots, cols, fov_id, condition_name, color, lw,
                               condition_output_dir, axis_limit)
        except Exception as e:
            print(f"    Error processing {fov_id}: {e}")
    print(f"  Completed: Spider plots saved to {condition_output_dir}")


def create_spider_plot(df_spots, cols, fov_id, condition_name, color, lw,
                        output_dir, axis_limit):
    fig, ax = plt.subplots(figsize=(10, 10))
    n_tracks = 0
    track_lengths = []

    for track_id, track_data in df_spots.groupby(cols['id']):
        if len(track_data) < 2:
            continue
        n_tracks += 1
        track_lengths.append(len(track_data))
        x0 = track_data[cols['x']].iloc[0]
        y0 = track_data[cols['y']].iloc[0]
        xn = track_data[cols['x']].values - x0
        yn = track_data[cols['y']].values - y0
        ax.plot(xn, yn, color=color, linewidth=lw, alpha=0.6, zorder=2)
        ax.plot(0, 0, 'o', color='green', markersize=4, alpha=0.3, zorder=1)
        ax.plot(xn[-1], yn[-1], 'o', color=color, markersize=3, alpha=0.5, zorder=3)

    _style_ax(ax, axis_limit)
    _add_reference_circles(ax, axis_limit)

    title = (f'{condition_name} - {fov_id} Spider Plot\n'
             f'N={n_tracks} tracks, Mean length={np.mean(track_lengths):.1f} points'
             if track_lengths else f'{condition_name} - {fov_id} (no tracks)')
    ax.set_title(title, fontsize=15, fontweight='bold', pad=15)

    stats_text = f'Tracks: {n_tracks}\n'
    if track_lengths:
        stats_text += f'Avg length: {np.mean(track_lengths):.1f} pts\n'
    stats_text += f'Axis: ±{axis_limit} µm'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8,
                      edgecolor=color, linewidth=2),
            fontsize=11, family='monospace')

    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='green',
               markersize=8, label='Start (origin)', alpha=0.7),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=color,
               markersize=6, label='End position', alpha=0.7),
        Line2D([0], [0], color=color, linewidth=lw, label='Track path', alpha=0.7),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=10, framealpha=0.9)
    plt.tight_layout()
    output_path = os.path.join(output_dir, f'{fov_id}_spider.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()


def create_combined_spider_plot(condition_path, condition_name, color, lw,
                                 output_dir, axis_limit=3.5):
    print(f"  Creating combined spider plot for all FOVs in {condition_name}...")
    spots_files = sorted(glob(os.path.join(condition_path, "FOV*-spots.csv")))
    if not spots_files:
        print("    WARNING: No spots files found")
        return

    fig, ax = plt.subplots(figsize=(12, 12))
    total_tracks = 0
    all_track_lengths = []
    all_end_distances = []

    for spots_file in spots_files:
        try:
            df_spots = pd.read_csv(spots_file)
            cols = detect_cols(df_spots)
            df_spots = df_spots.sort_values([cols['id'], cols['t']])
            for _, track_data in df_spots.groupby(cols['id']):
                if len(track_data) < 2:
                    continue
                total_tracks += 1
                all_track_lengths.append(len(track_data))
                x0 = track_data[cols['x']].iloc[0]
                y0 = track_data[cols['y']].iloc[0]
                xn = track_data[cols['x']].values - x0
                yn = track_data[cols['y']].values - y0
                all_end_distances.append(np.sqrt(xn[-1]**2 + yn[-1]**2))
                ax.plot(xn, yn, color=color, linewidth=lw * 0.67, alpha=0.3, zorder=2)
        except Exception as e:
            print(f"    Error: {e}")

    ax.plot(0, 0, 'o', color='green', markersize=8, alpha=0.8, zorder=5,
            label='Origin (all tracks)')
    _style_ax(ax, axis_limit)
    _add_reference_circles(ax, axis_limit, radii=[0.5, 1, 1.5, 2, 2.5, 3])

    title = (f'{condition_name} - Combined Spider Plot (All FOVs)\n'
             f'Total N={total_tracks} tracks from {len(spots_files)} FOVs')
    ax.set_title(title, fontsize=17, fontweight='bold', pad=20)

    stats_text = (f'Total tracks: {total_tracks}\n'
                  f'FOVs: {len(spots_files)}\n')
    if all_track_lengths:
        stats_text += (f'Avg track length: {np.mean(all_track_lengths):.1f} pts\n'
                       f'Mean end distance: {np.mean(all_end_distances):.2f} µm\n'
                       f'Median end distance: {np.median(all_end_distances):.2f} µm\n')
    stats_text += f'Axis range: ±{axis_limit} µm'
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.9,
                      edgecolor=color, linewidth=3),
            fontsize=12, family='monospace')
    ax.legend(loc='upper right', fontsize=12, framealpha=0.9)
    plt.tight_layout()
    output_path = os.path.join(output_dir, f'{condition_name}_COMBINED_spider.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"    Combined plot saved: {output_path}")


def create_all_conditions_comparison_spider(base_dir, conditions, output_dir, axis_limit=3.5):
    print(f"\n{'='*60}")
    print("Creating 2x2 comparison spider plot for all conditions")
    print(f"{'='*60}")

    fig, axes = plt.subplots(2, 2, figsize=(18, 18))
    fig.suptitle('Spider Plot Comparison - All Conditions',
                 fontsize=20, fontweight='bold', y=0.995)
    axes = axes.flatten()

    for idx, (condition_name, cfg) in enumerate(conditions.items()):
        ax = axes[idx]
        color = cfg['color']
        lw    = cfg['lw']
        condition_path = os.path.join(base_dir, condition_name)

        if not os.path.exists(condition_path):
            ax.text(0.5, 0.5, f'{condition_name}\nNo data found',
                    ha='center', va='center', fontsize=16)
            ax.set_title(condition_name, fontsize=15, fontweight='bold')
            continue

        spots_files = sorted(glob(os.path.join(condition_path, "FOV*-spots.csv")))
        if not spots_files:
            ax.text(0.5, 0.5, f'{condition_name}\nNo spots files',
                    ha='center', va='center', fontsize=16)
            ax.set_title(condition_name, fontsize=15, fontweight='bold')
            continue

        print(f"  Processing {condition_name}...")
        total_tracks = 0
        all_end_distances = []

        for spots_file in spots_files:
            try:
                df_spots = pd.read_csv(spots_file)
                cols = detect_cols(df_spots)
                df_spots = df_spots.sort_values([cols['id'], cols['t']])
                for _, track_data in df_spots.groupby(cols['id']):
                    if len(track_data) < 2:
                        continue
                    total_tracks += 1
                    x0 = track_data[cols['x']].iloc[0]
                    y0 = track_data[cols['y']].iloc[0]
                    xn = track_data[cols['x']].values - x0
                    yn = track_data[cols['y']].values - y0
                    all_end_distances.append(np.sqrt(xn[-1]**2 + yn[-1]**2))
                    ax.plot(xn, yn, color=color, linewidth=lw * 0.5,
                            alpha=0.25, zorder=2)
            except Exception:
                continue

        ax.plot(0, 0, 'o', color='green', markersize=6, alpha=0.8, zorder=5)
        _style_ax(ax, axis_limit)
        _add_reference_circles(ax, axis_limit)

        ax.set_title(condition_name, fontsize=16, fontweight='bold',
                     color=color, pad=10)
        stats_text = f'N = {total_tracks}\n'
        if all_end_distances:
            stats_text += f'Mean end: {np.mean(all_end_distances):.2f} µm'
        ax.text(0.02, 0.98, stats_text, transform=ax.transAxes,
                verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.85,
                          edgecolor=color, linewidth=2),
                fontsize=11, family='monospace')

    plt.tight_layout()
    output_path = os.path.join(output_dir, 'ALL_CONDITIONS_comparison_spider.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"\n  Comparison plot saved: {output_path}")


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main():
    print("\n" + "="*80)
    print("SPIDER PLOT GENERATOR v2 - Batch Processing for Multiple Conditions")
    print("="*80)

    base_dir   = input("\nEnter the base directory containing condition folders: ").strip().replace('"', '').replace("'", "")
    output_dir = input("Enter the output directory for spider plots: ").strip().replace('"', '').replace("'", "")

    axis_input = input("Enter axis limit in µm (default 3.5, press Enter for default): ").strip()
    axis_limit = float(axis_input) if axis_input else 3.5
    print(f"\nAxis range: ±{axis_limit} µm (total {axis_limit*2} x {axis_limit*2} µm)")

    n_top_input = input("Number of top tracks to show (default 10, press Enter for default): ").strip()
    n_top = int(n_top_input) if n_top_input else 10

    os.makedirs(output_dir, exist_ok=True)

    print("\n" + "="*80)
    print("PLOT OPTIONS:")
    print("  1. Individual spider plots for each FOV in each condition")
    print("  2. Combined spider plot for each condition (all FOVs merged)")
    print("  3. 2x2 comparison plot (one panel per condition)")
    print(f"  4. Top-{n_top} longest tracks by DURATION (one plot per track, per condition)")
    print(f"  5. Top-{n_top} longest tracks by PATH LENGTH (one plot per track, per condition)")
    print("  6. All of the above")
    print("="*80)

    plot_choice = input("\nEnter your choice (1-6, or comma-separated, e.g. 4,5): ").strip()
    choices = set(c.strip() for c in plot_choice.replace('-', ',').split(','))
    # '6' → all
    if '6' in choices:
        choices = {'1', '2', '3', '4', '5'}

    print("\n" + "="*80)
    print("PROCESSING CONDITIONS")
    print("="*80)

    for condition_name, cfg in CONDITIONS.items():
        color = cfg['color']
        lw    = cfg['lw']
        condition_path = os.path.join(base_dir, condition_name)

        if not os.path.exists(condition_path):
            print(f"\nWARNING: Directory not found: {condition_path}")
            continue

        if '1' in choices:
            generate_spider_plots_for_condition(condition_path, condition_name,
                                                color, lw, output_dir, axis_limit)
        if '2' in choices:
            create_combined_spider_plot(condition_path, condition_name,
                                        color, lw, output_dir, axis_limit)
        if '4' in choices:
            plot_top10_longest_duration(condition_path, condition_name,
                                        color, lw, output_dir, axis_limit, n_top)
        if '5' in choices:
            plot_top10_longest_path(condition_path, condition_name,
                                    color, lw, output_dir, axis_limit, n_top)

    if '3' in choices:
        create_all_conditions_comparison_spider(base_dir, CONDITIONS, output_dir, axis_limit)

    print("\n" + "="*80)
    print("SPIDER PLOT GENERATION COMPLETE!")
    print(f"Results saved to: {output_dir}")
    print("="*80)

    print("\nGenerated files:")
    if '1' in choices:
        print("  • Individual FOV spider plots in condition-specific subdirectories")
    if '2' in choices:
        print("  • Combined spider plots (one per condition)")
    if '3' in choices:
        print("  • 2x2 comparison plot with all conditions")
    if '4' in choices:
        print(f"  • Top-{n_top} duration tracks – individual spider plots per condition")
    if '5' in choices:
        print(f"  • Top-{n_top} path-length tracks – individual spider plots per condition")
    print("="*80)


if __name__ == "__main__":
    main()
