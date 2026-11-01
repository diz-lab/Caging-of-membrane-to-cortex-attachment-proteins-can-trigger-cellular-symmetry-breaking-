import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from glob import glob

def run_smt_analysis_v6():
    print("--- SMT Robust Analysis v6 (High-Diversity Dark Palette) ---")
    input_dir = input("Enter the input directory path: ").strip().replace('"', '').replace("'", "")
    output_dir = input("Enter the output directory path: ").strip().replace('"', '').replace("'", "")
    frame_interval = float(input("Enter the frame interval (dt) in seconds: "))
    spider_max_scale = float(input("Enter max scale for Spider Plot axes (microns): ") or 2)
    time_threshold = float(input("Enter minimum duration threshold (seconds): ") or 1.0)
    
    os.makedirs(output_dir, exist_ok=True)
    spot_files = sorted(glob(os.path.join(input_dir, "FOV*-spots.csv")))
    
    if not spot_files:
        print(f"ERROR: No files found in {input_dir}.")
        return

    # --- EXPANDED DARK/SATURATED PALETTE (No light colors) ---
    # Includes deep blues, greens, reds, purples, oranges, and teals
    dark_colors = [
        '#000080', '#8B0000', '#006400', '#2F4F4F', '#4B0082', '#191970', 
        '#800000', '#556B2F', '#0000CD', '#483D8B', '#008080', '#000000',
        '#800080', '#A52A2A', '#2E8B57', '#1E90FF', '#D2691E', '#C71585',
        '#4682B4', '#004d40', '#3e2723', '#1a237e', '#b71c1c', '#01579b'
    ]

    summary_data = []

    for s_file in spot_files:
        fov_id = os.path.basename(s_file).split('-')[0]
        print(f"\nProcessing {fov_id}...")
        
        try:
            df = pd.read_csv(s_file)
            
            # Detect Column Names
            cols = {
                'x': 'POSITION_X' if 'POSITION_X' in df.columns else 'x',
                'y': 'POSITION_Y' if 'POSITION_Y' in df.columns else 'y',
                't': 'POSITION_T' if 'POSITION_T' in df.columns else 't',
                'id': 'TRACK_ID' if 'TRACK_ID' in df.columns else 'track_id'
            }
            
            df = df.sort_values([cols['id'], cols['t']])
            
            # 1. Temporal Filter
            track_stats = df.groupby(cols['id']).agg(
                Duration=(cols['t'], lambda x: x.max() - x.min())
            )
            initial_count = len(track_stats)
            valid_ids = track_stats[track_stats['Duration'] >= time_threshold].index
            df_filtered = df[df[cols['id']].isin(valid_ids)].copy()
            
            if df_filtered.empty:
                print(f"Skipping {fov_id}: No tracks meet time threshold.")
                continue

            # 2. Dynamics Calculation
            def process_track(group):
                dr = np.sqrt(group[cols['x']].diff()**2 + group[cols['y']].diff()**2)
                t_rel = group[cols['t']] - group[cols['t']].iloc[0]
                return pd.DataFrame({'Disp': dr, 'T_Rel': t_rel, 'TRACK_ID': group[cols['id']]})

            dynamics = df_filtered.groupby(cols['id'], group_keys=False).apply(process_track)
            df_filtered['Disp'] = dynamics['Disp']
            df_filtered['T_Rel'] = dynamics['T_Rel']
            
            # 3. Distribution Plot
            track_medians = df_filtered.groupby(cols['id'])['Disp'].median().dropna()
            clean_medians = track_medians[np.isfinite(track_medians) & (track_medians < 50)]

            plt.figure(figsize=(8, 8))
            if not clean_medians.empty:
                sns.histplot(clean_medians, bins=100, kde=True, color='gray', alpha=0.3)
                plt.axvline(clean_medians.mean(), color='blue', ls='--', label=f'Mean: {clean_medians.mean():.3f}')
                plt.axvline(clean_medians.median(), color='green', ls='-', label=f'Med: {clean_medians.median():.3f}')
                plt.axvline(0.2, color='red', ls=':', label='200nm Threshold')
                plt.title(f"{fov_id} Displacement Distribution")
                plt.xlabel("Median Displacement (µm)")
                plt.legend()
            plt.savefig(os.path.join(output_dir, f"{fov_id}_dist.png"), dpi=300)
            plt.close()

            # 4. Spatial & Spider Plots (All)
            plot_spatial(df_filtered, f"{fov_id}_All", output_dir, dark_colors, cols)
            plot_spider(df_filtered, f"{fov_id}_All", output_dir, spider_max_scale, dark_colors, cols)

            # 5. Mobile Filter (>= 200nm)
            mobile_ids = track_medians[track_medians >= 0.2].index
            df_mobile = df_filtered[df_filtered[cols['id']].isin(mobile_ids)].copy()
            
            if not df_mobile.empty:
                # 6. Re-plot Mobile
                plot_spatial(df_mobile, f"{fov_id}_Mobile", output_dir, dark_colors, cols)
                plot_spider(df_mobile, f"{fov_id}_Mobile", output_dir, spider_max_scale, dark_colors, cols)
                
                # 7. Disp vs Time Plot
                plt.figure(figsize=(8, 8))
                for i, (tid, group) in enumerate(df_mobile.groupby(cols['id'])):
                    c = dark_colors[i % len(dark_colors)]
                    v = group.dropna(subset=['Disp'])
                    plt.plot(v['T_Rel'], v['Disp'], '-o', color=c, ms=3, alpha=0.6, lw=0.8)
                plt.title(f"{fov_id} Mobile: Displacement vs Time Interval")
                plt.xlabel("Time (s)"); plt.ylabel("Inst. Displacement (µm)")
                plt.xlim(0, df_mobile['T_Rel'].max())
                plt.savefig(os.path.join(output_dir, f"{fov_id}_mobile_dynamics.png"), dpi=300)
                plt.close()

            summary_data.append({
                'FOV': fov_id, 'Initial_Tracks': initial_count, 
                'Below_200nm': (len(valid_ids) - len(mobile_ids)), 'Remaining_Mobile': len(mobile_ids)
            })

        except Exception as e:
            print(f"Error on {fov_id}: {e}")

    pd.DataFrame(summary_data).to_csv(os.path.join(output_dir, "report.csv"), index=False)
    print("\nProcessing Complete.")

# --- Plotting Functions ---

def plot_spatial(df, label, out, colors, cols):
    plt.figure(figsize=(8, 8))
    for i, (tid, group) in enumerate(df.groupby(cols['id'])):
        c = colors[i % len(colors)]
        plt.plot(group[cols['x']], group[cols['y']], color=c, lw=1.0, alpha=0.8)
    plt.title(f"{label} Spatial Paths")
    plt.xlabel("X (µm)"); plt.ylabel("Y (µm)")
    plt.xlim(df[cols['x']].min(), df[cols['x']].max())
    plt.ylim(df[cols['y']].min(), df[cols['y']].max())
    plt.gca().set_aspect('equal')
    plt.grid(True, linestyle=':', alpha=0.3)
    plt.savefig(os.path.join(out, f"{label}_spatial.png"), dpi=300)
    plt.close()

def plot_spider(df, label, out, scale, colors, cols):
    plt.figure(figsize=(8, 8))
    for i, (tid, group) in enumerate(df.groupby(cols['id'])):
        c = colors[i % len(colors)]
        x0 = group[cols['x']] - group[cols['x']].iloc[0]
        y0 = group[cols['y']] - group[cols['y']].iloc[0]
        plt.plot(x0, y0, color=c, lw=1.0, alpha=0.8)
    plt.title(f"{label} Spider Plot")
    plt.xlabel("dX (µm)"); plt.ylabel("dY (µm)")
    plt.xlim(-scale, scale); plt.ylim(-scale, scale)
    plt.axhline(0, c='black', lw=0.5, alpha=0.3)
    plt.axvline(0, c='black', lw=0.5, alpha=0.3)
    plt.gca().set_aspect('equal')
    plt.grid(True, linestyle=':', alpha=0.3)
    plt.savefig(os.path.join(out, f"{label}_spider.png"), dpi=300)
    plt.close()

if __name__ == "__main__":
    run_smt_analysis_v6()
