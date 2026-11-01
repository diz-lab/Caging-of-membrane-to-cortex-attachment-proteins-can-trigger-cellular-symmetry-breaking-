import pandas as pd
import numpy as np
import os
from glob import glob

# =============================================================================
# CONFIGURATION
# =============================================================================

FOLDER_FRAME_RATES = {
    "MoeL":       0.02,
    "MoeLActin":  0.02,
    "MoeS":       0.03,
    "MoeSActin":  0.03,
}

DISPLACEMENT_THRESHOLD = 0.2   # µm

# =============================================================================
# HELPERS
# =============================================================================

def detect_cols(df):
    return {
        'x':  'POSITION_X' if 'POSITION_X' in df.columns else 'x',
        'y':  'POSITION_Y' if 'POSITION_Y' in df.columns else 'y',
        't':  'POSITION_T' if 'POSITION_T' in df.columns else 't',
        'id': 'TRACK_ID'   if 'TRACK_ID'   in df.columns else 'track_id',
    }

def fmt(v):
    return f"{v:.3f}".rstrip('0').rstrip('.')

def save_pair(spots_df, tracks_df, track_ids, spots_col_id, tracks_col_id,
              folder_path, filename_stem):
    sp = spots_df[spots_df[spots_col_id].isin(track_ids)].copy()
    tr = tracks_df[tracks_df[tracks_col_id].isin(track_ids)].copy()
    sp.to_csv(os.path.join(folder_path, f"{filename_stem}-mobile spots.csv"),  index=False)
    tr.to_csv(os.path.join(folder_path, f"{filename_stem}-mobile tracks.csv"), index=False)

# =============================================================================
# PROCESS ONE FOV
# =============================================================================

def process_fov(folder_path, fov_id, tracks_df, spots_df, comprehensive_df,
                dur_thresholds, frame_interval):
    """
    Pipeline per FOV:
      1. From the comprehensive results file, keep only tracks labelled 'Mobile'
      2. Filter those tracks to median displacement >= 0.2 µm
      3. Filter that pool by each duration threshold  →  save tracks + spots
      4. Return a summary row dict
    """
    cols       = detect_cols(spots_df)
    track_cols = detect_cols(tracks_df)
    folder_name = os.path.basename(folder_path)

    spots_df = spots_df.sort_values([cols['id'], cols['t']])
    n_total = spots_df[cols['id']].nunique()

    # ------------------------------------------------------------------
    # 1. Classification filter: keep only 'Mobile' tracks from
    #    the comprehensive results file
    # ------------------------------------------------------------------
    if comprehensive_df is None:
        print(f"    [WARNING] No comprehensive results file found – using all tracks.")
        classified_ids = spots_df[cols['id']].unique()
        n_classified = len(classified_ids)
    else:
        # Detect the classification column (case-insensitive)
        class_col = next(
            (c for c in comprehensive_df.columns if c.strip().lower() == 'classification'),
            None
        )
        id_col_comp = next(
            (c for c in comprehensive_df.columns if c.strip().upper() == 'TRACK_ID'),
            None
        )
        if class_col is None or id_col_comp is None:
            print(f"    [WARNING] Could not find TRACK_ID or classification column in "
                  f"comprehensive file. Columns: {list(comprehensive_df.columns)}")
            classified_ids = spots_df[cols['id']].unique()
            n_classified = len(classified_ids)
        else:
            mobile_class = comprehensive_df[
                comprehensive_df[class_col].str.strip().str.lower() == 'mobile'
            ][id_col_comp]
            classified_ids = mobile_class.values
            n_classified = len(classified_ids)

    print(f"    Total tracks in file          : {n_total}")
    print(f"    Labelled 'Mobile' in comp file: {n_classified}")

    # Restrict spots and tracks to classified-mobile IDs
    spots_class  = spots_df[spots_df[cols['id']].isin(classified_ids)].copy()
    tracks_class = tracks_df[tracks_df[track_cols['id']].isin(classified_ids)].copy()

    if spots_class.empty:
        print(f"    [INFO] No 'Mobile' tracks found – skipping.")
        return None

    # ------------------------------------------------------------------
    # 2. Per-step displacement → per-track median on the classified pool
    # ------------------------------------------------------------------
    def step_disp(group):
        return np.sqrt(group[cols['x']].diff()**2 + group[cols['y']].diff()**2)

    spots_class = spots_class.copy()
    spots_class['step_disp'] = (
        spots_class.groupby(cols['id'], group_keys=False).apply(step_disp)
    )
    track_medians = spots_class.groupby(cols['id'])['step_disp'].median().dropna()

    print(f"    Displacement stats on Mobile pool ({len(track_medians)} tracks):")
    print(f"      min={track_medians.min():.4f}  median={track_medians.median():.4f}"
          f"  mean={track_medians.mean():.4f}  max={track_medians.max():.4f}")
    below = (track_medians < DISPLACEMENT_THRESHOLD).sum()
    print(f"      Below {fmt(DISPLACEMENT_THRESHOLD)} µm: {below}  |  "
          f"Above: {len(track_medians) - below}")

    mobile_ids = track_medians[track_medians >= DISPLACEMENT_THRESHOLD].index
    n_mobile   = len(mobile_ids)

    summary = {
        'Folder':                                       folder_name,
        'FOV':                                          fov_id,
        'Total tracks':                                 n_total,
        "Labelled 'Mobile' (comp file)":                n_classified,
        f'Disp >= {fmt(DISPLACEMENT_THRESHOLD)} µm':   n_mobile,
    }

    spots_clean = spots_class.drop(columns=['step_disp'])

    if n_mobile == 0:
        print(f"    Disp >= {fmt(DISPLACEMENT_THRESHOLD)} µm  →  0 tracks")
        for dt in dur_thresholds:
            summary[f'Disp >= {fmt(DISPLACEMENT_THRESHOLD)} µm + dur >= {fmt(dt)}s'] = 0
        return summary

    save_pair(spots_clean, tracks_class, mobile_ids,
              cols['id'], track_cols['id'],
              folder_path, f"{fov_id}-{fmt(DISPLACEMENT_THRESHOLD)}um")
    print(f"    Disp >= {fmt(DISPLACEMENT_THRESHOLD)} µm  →  {n_mobile} tracks"
          f"  |  {fov_id}-{fmt(DISPLACEMENT_THRESHOLD)}um-mobile tracks/spots.csv")

    # ------------------------------------------------------------------
    # 3. Duration sort on the displacement-filtered mobile pool
    # ------------------------------------------------------------------

    # Prefer TRACK_DURATION from tracks file; fall back to computing from POSITION_T
    if 'TRACK_DURATION' in tracks_class.columns:
        track_durations = tracks_class.set_index(track_cols['id'])['TRACK_DURATION']
        dur_source = "TRACK_DURATION column"
    else:
        print(f"    [NOTE] TRACK_DURATION not in tracks file – computing from POSITION_T.")
        raw_dur = spots_class.groupby(cols['id'])[cols['t']].agg(lambda x: x.max() - x.min())
        median_raw = raw_dur.median()
        if median_raw > 500:
            track_durations = (raw_dur * frame_interval).rename('Duration')
            print(f"    [NOTE] Frame numbers detected (median={median_raw:.1f}); "
                  f"converting with dt={frame_interval} s.")
        else:
            track_durations = raw_dur.rename('Duration')
            print(f"    [NOTE] POSITION_T in seconds (median={median_raw:.4f} s).")
        dur_source = "POSITION_T (derived)"

    mob_durations = track_durations[track_durations.index.isin(mobile_ids)]
    print(f"    Duration stats on {n_mobile} disp-filtered tracks ({dur_source}):")
    print(f"      min={mob_durations.min():.4f}  median={mob_durations.median():.4f}"
          f"  max={mob_durations.max():.4f}")

    print(f"    --- Duration sort ---")
    for dt in dur_thresholds:
        dur_ids = mob_durations[mob_durations >= dt].index
        n_dur   = len(dur_ids)
        col_name = f'Disp >= {fmt(DISPLACEMENT_THRESHOLD)} µm + dur >= {fmt(dt)}s'
        summary[col_name] = n_dur

        if n_dur > 0:
            save_pair(spots_clean, tracks_class, dur_ids,
                      cols['id'], track_cols['id'],
                      folder_path, f"{fov_id}-{fmt(DISPLACEMENT_THRESHOLD)}um-dur{fmt(dt)}s")
            print(f"      dur >= {fmt(dt)}s  →  {n_dur} tracks"
                  f"  |  {fov_id}-{fmt(DISPLACEMENT_THRESHOLD)}um-dur{fmt(dt)}s-mobile tracks/spots.csv")
        else:
            print(f"      dur >= {fmt(dt)}s  →  0 tracks")

    return summary

# =============================================================================
# PROCESS ONE FOLDER
# =============================================================================

def process_folder(folder_path, frame_interval, dur_thresholds):
    folder_name = os.path.basename(folder_path)
    spot_files  = sorted(glob(os.path.join(folder_path, "FOV*-spots.csv")))

    if not spot_files:
        print(f"  [WARNING] No FOV*-spots.csv files found in: {folder_path}")
        return []

    summaries = []
    for s_file in spot_files:
        fov_id = os.path.basename(s_file).split('-')[0]
        print(f"\n  {folder_name} / {fov_id}")
        print(f"  {'─' * 40}")

        t_file = os.path.join(folder_path, f"{fov_id}-tracks.csv")
        if not os.path.exists(t_file):
            print(f"    [WARNING] Tracks file not found – skipping.")
            continue

        # Look for the comprehensive results file for this FOV
        comp_pattern = os.path.join(folder_path, f"{fov_id}_comprehensive_results.csv")
        comp_files   = glob(comp_pattern)
        if comp_files:
            comprehensive_df = pd.read_csv(comp_files[0])
            print(f"    Comprehensive file: {os.path.basename(comp_files[0])}")
        else:
            comprehensive_df = None
            print(f"    [NOTE] No comprehensive results file found for {fov_id}.")

        try:
            spots_df  = pd.read_csv(s_file)
            tracks_df = pd.read_csv(t_file)
            row = process_fov(folder_path, fov_id, tracks_df, spots_df,
                              comprehensive_df, dur_thresholds, frame_interval)
            if row:
                summaries.append(row)
        except Exception as e:
            print(f"    [ERROR] {e}")

    return summaries

# =============================================================================
# PROMPTS
# =============================================================================

def prompt_duration_thresholds():
    while True:
        raw = input("How many duration thresholds would you like? (1, 2, or 3): ").strip()
        if raw in ('1', '2', '3'):
            n = int(raw)
            break
        print("  Please enter 1, 2, or 3.")

    thresholds = []
    for i in range(1, n + 1):
        while True:
            raw = input(f"  Duration threshold {i} (seconds): ").strip()
            try:
                thresholds.append(float(raw))
                break
            except ValueError:
                print("  Please enter a valid number.")

    return sorted(set(thresholds))

# =============================================================================
# MAIN
# =============================================================================

def run_batch():
    base_dir = input("Enter the base directory containing the 4 condition folders: ").strip().strip('"').strip("'")

    print()
    dur_thresholds = prompt_duration_thresholds()

    print()
    print("=" * 70)
    print("  Batch Mobile Track Filter")
    print(f"  Base directory          : {base_dir}")
    print(f"  Step 1 – Classification : 'Mobile' label from FOV*_comprehensive_results.csv")
    print(f"  Step 2 – Displacement   : >= {fmt(DISPLACEMENT_THRESHOLD)} µm")
    print(f"  Step 3 – Duration       : {dur_thresholds} s")
    print("=" * 70)

    all_summaries = []

    for folder_name, frame_interval in FOLDER_FRAME_RATES.items():
        folder_path = os.path.join(base_dir, folder_name)

        if not os.path.isdir(folder_path):
            print(f"\n[WARNING] Folder not found, skipping: {folder_path}")
            continue

        print(f"\n{'═' * 70}")
        print(f"  Folder : {folder_name}  (dt = {frame_interval} s)")
        print(f"{'═' * 70}")

        rows = process_folder(folder_path, frame_interval, dur_thresholds)
        all_summaries.extend(rows)

    # ------------------------------------------------------------------
    # Save summary tables
    # ------------------------------------------------------------------
    if all_summaries:
        summary_df = pd.DataFrame(all_summaries)

        for folder_name in FOLDER_FRAME_RATES:
            sub = summary_df[summary_df['Folder'] == folder_name]
            if sub.empty:
                continue
            out = os.path.join(base_dir, folder_name, f"{folder_name}-summary.csv")
            sub.to_csv(out, index=False)
            print(f"\n  Folder summary saved : {out}")

        combined = os.path.join(base_dir, "all_conditions_summary.csv")
        summary_df.to_csv(combined, index=False)

        print()
        print("=" * 70)
        print("  SUMMARY TABLE")
        print("=" * 70)
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', 200)
        print(summary_df.to_string(index=False))
        print()
        print(f"  Combined summary saved : {combined}")

    print()
    print("=" * 70)
    print("  Batch processing complete.")
    print("=" * 70)


if __name__ == "__main__":
    run_batch()
