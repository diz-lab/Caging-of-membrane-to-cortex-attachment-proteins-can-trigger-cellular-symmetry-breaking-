import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from glob import glob
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde
import warnings
warnings.filterwarnings('ignore')

def gaussian(x, amplitude, mean, stddev):
    """Gaussian function for fitting"""
    return amplitude * np.exp(-((x - mean) ** 2) / (2 * stddev ** 2))

def calculate_r_squared(y_actual, y_predicted):
    """Calculate R-squared value"""
    ss_res = np.sum((y_actual - y_predicted) ** 2)
    ss_tot = np.sum((y_actual - np.mean(y_actual)) ** 2)
    return 1 - (ss_res / ss_tot) if ss_tot != 0 else 0

def calculate_first_passage_time(track_data, cols, target_radius):
    """
    Calculate first passage time: time to first exit from initial position
    within a given radius
    """
    x0 = track_data[cols['x']].iloc[0]
    y0 = track_data[cols['y']].iloc[0]
    t0 = track_data[cols['t']].iloc[0]
    
    # Calculate distance from origin for each point
    distances = np.sqrt((track_data[cols['x']] - x0)**2 + (track_data[cols['y']] - y0)**2)
    
    # Find first time the particle exits the radius
    exit_points = distances > target_radius
    if exit_points.any():
        first_exit_idx = exit_points.idxmax()
        exit_time = track_data.loc[first_exit_idx, cols['t']] - t0
        return exit_time
    else:
        return np.nan

def calculate_free_exit_time(track_data, cols, confinement_threshold=0.5):
    """
    Calculate free exit time: duration of confined motion before transitioning to free diffusion
    Using displacement-based criterion
    """
    # Calculate instantaneous displacements
    displacements = np.sqrt(track_data[cols['x']].diff()**2 + track_data[cols['y']].diff()**2)
    
    # Find transition point where displacement consistently exceeds threshold
    window_size = 3
    if len(displacements) < window_size:
        return np.nan
    
    rolling_mean = displacements.rolling(window=window_size).mean()
    transition_idx = (rolling_mean > confinement_threshold).idxmax()
    
    if pd.isna(transition_idx) or transition_idx == displacements.index[0]:
        # No clear transition found, return full track duration
        return track_data[cols['t']].iloc[-1] - track_data[cols['t']].iloc[0]
    
    exit_time = track_data.loc[transition_idx, cols['t']] - track_data[cols['t']].iloc[0]
    return exit_time if exit_time > 0 else np.nan

def process_condition(condition_path, condition_name, colors_dict, radii=[1, 2, 3]):
    """
    Process all FOV files for a single condition
    """
    print(f"\n{'='*60}")
    print(f"Processing condition: {condition_name}")
    print(f"{'='*60}")
    
    # Find all spots and tracks files
    spots_files = sorted(glob(os.path.join(condition_path, "FOV*-spots.csv")))
    tracks_files = sorted(glob(os.path.join(condition_path, "FOV*-tracks.csv")))
    
    if not spots_files:
        print(f"WARNING: No spots files found in {condition_path}")
        return None
    
    print(f"Found {len(spots_files)} FOV files")
    
    # Storage for all metrics
    all_fpt_data = {radius: [] for radius in radii}
    all_lifetimes = []
    all_end_to_end = []
    all_path_lengths = []
    all_free_exit_times = []
    
    # Process each FOV
    for spots_file in spots_files:
        fov_id = os.path.basename(spots_file).split('-')[0]
        print(f"  Processing {fov_id}...")
        
        try:
            # Load spots data
            df_spots = pd.read_csv(spots_file)
            
            # Detect column names
            cols = {
                'x': 'POSITION_X' if 'POSITION_X' in df_spots.columns else 'x',
                'y': 'POSITION_Y' if 'POSITION_Y' in df_spots.columns else 'y',
                't': 'POSITION_T' if 'POSITION_T' in df_spots.columns else 't',
                'id': 'TRACK_ID' if 'TRACK_ID' in df_spots.columns else 'track_id'
            }
            
            # Sort by track and time
            df_spots = df_spots.sort_values([cols['id'], cols['t']])

            # --- FIX: Clean time column to prevent negative FPT values ---
            # Reset time to 0 for each track so t0 is always 0
            def reset_track_time(group):
                group = group.copy()
                group[cols['t']] = group[cols['t']] - group[cols['t']].iloc[0]
                return group

            df_spots = df_spots.groupby(cols['id'], group_keys=False).apply(reset_track_time)

            # Remove any tracks where time is still non-monotonic after sorting
            def time_is_monotonic(group):
                times = group[cols['t']].values
                return all(times[i] <= times[i+1] for i in range(len(times) - 1))

            df_spots = df_spots.groupby(cols['id'], group_keys=False).filter(time_is_monotonic)
            # --- END FIX ---

            # Process each track
            for track_id, track_data in df_spots.groupby(cols['id']):
                if len(track_data) < 3:  # Skip very short tracks
                    continue
                
                # 1. First Passage Time for different radii
                for radius in radii:
                    fpt = calculate_first_passage_time(track_data, cols, radius)
                    if not np.isnan(fpt) and fpt > 0:
                        all_fpt_data[radius].append(fpt)
                
                # 2. Track lifetime
                lifetime = track_data[cols['t']].iloc[-1] - track_data[cols['t']].iloc[0]
                all_lifetimes.append(lifetime)
                
                # 3. End-to-end distance
                x0, y0 = track_data[cols['x']].iloc[0], track_data[cols['y']].iloc[0]
                xf, yf = track_data[cols['x']].iloc[-1], track_data[cols['y']].iloc[-1]
                end_to_end = np.sqrt((xf - x0)**2 + (yf - y0)**2)
                all_end_to_end.append(end_to_end)
                
                # 4. Path length (total distance traveled)
                dx = track_data[cols['x']].diff()
                dy = track_data[cols['y']].diff()
                path_length = np.sqrt(dx**2 + dy**2).sum()
                all_path_lengths.append(path_length)
                
                # 5. Free Exit Time
                fet = calculate_free_exit_time(track_data, cols)
                if not np.isnan(fet) and fet > 0:
                    all_free_exit_times.append(fet)
        
        except Exception as e:
            print(f"    Error processing {fov_id}: {e}")
            continue
    
    # Package results
    results = {
        'condition_name': condition_name,
        'color': colors_dict[condition_name],
        'fpt_data': all_fpt_data,
        'lifetimes': all_lifetimes,
        'end_to_end': all_end_to_end,
        'path_lengths': all_path_lengths,
        'free_exit_times': all_free_exit_times,
        'n_tracks': len(all_lifetimes)
    }
    
    print(f"  Completed: {len(all_lifetimes)} tracks analyzed")
    
    return results

def plot_fpt_distributions(all_results, radii, output_dir):
    """
    Plot First Passage Time distributions for each radius
    """
    print("\nGenerating First Passage Time distributions...")
    
    for radius in radii:
        fig, axes = plt.subplots(2, 2, figsize=(14, 12))
        fig.suptitle(f'First Passage Time Distribution (Radius = {radius} µm)', 
                     fontsize=16, fontweight='bold')
        
        axes = axes.flatten()
        stats_summary = []
        
        for idx, result in enumerate(all_results):
            if result is None:
                continue
            
            ax = axes[idx]
            data = result['fpt_data'][radius]
            
            if len(data) == 0:
                ax.text(0.5, 0.5, 'No Data', ha='center', va='center', fontsize=14)
                ax.set_title(result['condition_name'])
                continue
            
            data = np.array(data)
            data = data[np.isfinite(data)]
            
            # Plot histogram
            counts, bins, patches = ax.hist(data, bins=50, alpha=0.6, 
                                           color=result['color'], 
                                           edgecolor='black', linewidth=0.5,
                                           density=True, label='Data')
            
            # Fit Gaussian
            try:
                bin_centers = (bins[:-1] + bins[1:]) / 2
                
                # Initial parameter guesses
                amplitude_guess = counts.max()
                mean_guess = np.mean(data)
                stddev_guess = np.std(data)
                
                popt, _ = curve_fit(gaussian, bin_centers, counts, 
                                   p0=[amplitude_guess, mean_guess, stddev_guess],
                                   maxfev=10000)
                
                # Generate fitted curve
                x_fit = np.linspace(data.min(), data.max(), 200)
                y_fit = gaussian(x_fit, *popt)
                
                ax.plot(x_fit, y_fit, 'r-', linewidth=3, label='Gaussian Fit')
                
                # Calculate R²
                y_hist_fit = gaussian(bin_centers, *popt)
                r_squared = calculate_r_squared(counts, y_hist_fit)
                
                # Extract parameters
                mean_fit = popt[1]
                std_fit = abs(popt[2])
                
                # Add statistics to plot
                stats_text = f'μ = {mean_fit:.3f} s\nσ = {std_fit:.3f} s\nR² = {r_squared:.3f}\nn = {len(data)}'
                ax.text(0.65, 0.95, stats_text, transform=ax.transAxes,
                       verticalalignment='top', bbox=dict(boxstyle='round', 
                       facecolor='wheat', alpha=0.8), fontsize=10)
                
                stats_summary.append({
                    'Condition': result['condition_name'],
                    'Radius_um': radius,
                    'Mean_s': mean_fit,
                    'StdDev_s': std_fit,
                    'R_squared': r_squared,
                    'N_tracks': len(data)
                })
                
            except Exception as e:
                print(f"    Could not fit Gaussian for {result['condition_name']}: {e}")
                mean_data = np.mean(data)
                std_data = np.std(data)
                stats_text = f'μ = {mean_data:.3f} s\nσ = {std_data:.3f} s\nn = {len(data)}'
                ax.text(0.65, 0.95, stats_text, transform=ax.transAxes,
                       verticalalignment='top', bbox=dict(boxstyle='round', 
                       facecolor='wheat', alpha=0.8), fontsize=10)
            
            ax.set_xlabel('First Passage Time (s)', fontsize=11)
            ax.set_ylabel('Probability Density', fontsize=11)
            ax.set_title(result['condition_name'], fontsize=13, fontweight='bold')
            ax.legend(loc='upper right')
            ax.grid(True, alpha=0.3, linestyle='--')
        
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, f'FPT_distribution_radius_{radius}um.png'), 
                   dpi=300, bbox_inches='tight')
        plt.close()
        
        # Save statistics
        if stats_summary:
            pd.DataFrame(stats_summary).to_csv(
                os.path.join(output_dir, f'FPT_statistics_radius_{radius}um.csv'), 
                index=False)

def plot_lifetime_distributions(all_results, output_dir):
    """
    Plot track lifetime distributions
    """
    print("\nGenerating lifetime distributions...")
    
    # Determine global bin range for consistency
    all_lifetimes = []
    for result in all_results:
        if result is not None:
            all_lifetimes.extend(result['lifetimes'])
    
    if not all_lifetimes:
        print("No lifetime data available")
        return
    
    bin_min = 0
    bin_max = np.percentile(all_lifetimes, 99)
    bins = np.linspace(bin_min, bin_max, 50)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle('Track Lifetime Distribution', fontsize=16, fontweight='bold')
    
    axes = axes.flatten()
    
    for idx, result in enumerate(all_results):
        if result is None:
            continue
        
        ax = axes[idx]
        data = np.array(result['lifetimes'])
        
        ax.hist(data, bins=bins, alpha=0.7, color=result['color'], 
               edgecolor='black', linewidth=0.5)
        
        mean_val = np.mean(data)
        median_val = np.median(data)
        
        ax.axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_val:.2f}s')
        ax.axvline(median_val, color='blue', linestyle='--', linewidth=2, label=f'Median: {median_val:.2f}s')
        
        stats_text = f'n = {len(data)}\nμ = {mean_val:.2f} s\nσ = {np.std(data):.2f} s'
        ax.text(0.65, 0.95, stats_text, transform=ax.transAxes,
               verticalalignment='top', bbox=dict(boxstyle='round', 
               facecolor='wheat', alpha=0.8), fontsize=10)
        
        ax.set_xlabel('Lifetime (s)', fontsize=11)
        ax.set_ylabel('Count', fontsize=11)
        ax.set_title(result['condition_name'], fontsize=13, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'lifetime_distributions.png'), 
               dpi=300, bbox_inches='tight')
    plt.close()

def plot_end_to_end_distributions(all_results, output_dir):
    """
    Plot end-to-end distance distributions
    """
    print("\nGenerating end-to-end distance distributions...")
    
    # Determine global bin range
    all_distances = []
    for result in all_results:
        if result is not None:
            all_distances.extend(result['end_to_end'])
    
    if not all_distances:
        print("No end-to-end distance data available")
        return
    
    bin_max = np.percentile(all_distances, 99)
    bins = np.linspace(0, bin_max, 50)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle('End-to-End Distance Distribution', fontsize=16, fontweight='bold')
    
    axes = axes.flatten()
    
    for idx, result in enumerate(all_results):
        if result is None:
            continue
        
        ax = axes[idx]
        data = np.array(result['end_to_end'])
        
        ax.hist(data, bins=bins, alpha=0.7, color=result['color'], 
               edgecolor='black', linewidth=0.5)
        
        mean_val = np.mean(data)
        median_val = np.median(data)
        
        ax.axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_val:.2f}µm')
        ax.axvline(median_val, color='blue', linestyle='--', linewidth=2, label=f'Median: {median_val:.2f}µm')
        
        stats_text = f'n = {len(data)}\nμ = {mean_val:.2f} µm\nσ = {np.std(data):.2f} µm'
        ax.text(0.65, 0.95, stats_text, transform=ax.transAxes,
               verticalalignment='top', bbox=dict(boxstyle='round', 
               facecolor='wheat', alpha=0.8), fontsize=10)
        
        ax.set_xlabel('End-to-End Distance (µm)', fontsize=11)
        ax.set_ylabel('Count', fontsize=11)
        ax.set_title(result['condition_name'], fontsize=13, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'end_to_end_distributions.png'), 
               dpi=300, bbox_inches='tight')
    plt.close()

def plot_path_length_distributions(all_results, output_dir):
    """
    Plot path length distributions
    """
    print("\nGenerating path length distributions...")
    
    # Determine global bin range
    all_lengths = []
    for result in all_results:
        if result is not None:
            all_lengths.extend(result['path_lengths'])
    
    if not all_lengths:
        print("No path length data available")
        return
    
    bin_max = np.percentile(all_lengths, 99)
    bins = np.linspace(0, bin_max, 50)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle('Path Length Distribution', fontsize=16, fontweight='bold')
    
    axes = axes.flatten()
    
    for idx, result in enumerate(all_results):
        if result is None:
            continue
        
        ax = axes[idx]
        data = np.array(result['path_lengths'])
        
        ax.hist(data, bins=bins, alpha=0.7, color=result['color'], 
               edgecolor='black', linewidth=0.5)
        
        mean_val = np.mean(data)
        median_val = np.median(data)
        
        ax.axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_val:.2f}µm')
        ax.axvline(median_val, color='blue', linestyle='--', linewidth=2, label=f'Median: {median_val:.2f}µm')
        
        stats_text = f'n = {len(data)}\nμ = {mean_val:.2f} µm\nσ = {np.std(data):.2f} µm'
        ax.text(0.65, 0.95, stats_text, transform=ax.transAxes,
               verticalalignment='top', bbox=dict(boxstyle='round', 
               facecolor='wheat', alpha=0.8), fontsize=10)
        
        ax.set_xlabel('Path Length (µm)', fontsize=11)
        ax.set_ylabel('Count', fontsize=11)
        ax.set_title(result['condition_name'], fontsize=13, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'path_length_distributions.png'), 
               dpi=300, bbox_inches='tight')
    plt.close()

def plot_free_exit_time_distributions(all_results, output_dir):
    """
    Plot Free Exit Time distributions
    """
    print("\nGenerating Free Exit Time distributions...")
    
    # Determine global bin range
    all_fet = []
    for result in all_results:
        if result is not None:
            all_fet.extend(result['free_exit_times'])
    
    if not all_fet:
        print("No Free Exit Time data available")
        return
    
    bin_max = np.percentile(all_fet, 99)
    bins = np.linspace(0, bin_max, 50)
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle('Free Exit Time Distribution', fontsize=16, fontweight='bold')
    
    axes = axes.flatten()
    
    for idx, result in enumerate(all_results):
        if result is None:
            continue
        
        ax = axes[idx]
        data = np.array(result['free_exit_times'])
        
        if len(data) == 0:
            ax.text(0.5, 0.5, 'No Data', ha='center', va='center', fontsize=14)
            ax.set_title(result['condition_name'])
            continue
        
        ax.hist(data, bins=bins, alpha=0.7, color=result['color'], 
               edgecolor='black', linewidth=0.5)
        
        mean_val = np.mean(data)
        median_val = np.median(data)
        
        ax.axvline(mean_val, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_val:.2f}s')
        ax.axvline(median_val, color='blue', linestyle='--', linewidth=2, label=f'Median: {median_val:.2f}s')
        
        stats_text = f'n = {len(data)}\nμ = {mean_val:.2f} s\nσ = {np.std(data):.2f} s'
        ax.text(0.65, 0.95, stats_text, transform=ax.transAxes,
               verticalalignment='top', bbox=dict(boxstyle='round', 
               facecolor='wheat', alpha=0.8), fontsize=10)
        
        ax.set_xlabel('Free Exit Time (s)', fontsize=11)
        ax.set_ylabel('Count', fontsize=11)
        ax.set_title(result['condition_name'], fontsize=13, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'free_exit_time_distributions.png'), 
               dpi=300, bbox_inches='tight')
    plt.close()

def generate_summary_report(all_results, output_dir):
    """
    Generate a comprehensive summary report
    """
    print("\nGenerating summary report...")
    
    summary_data = []
    
    for result in all_results:
        if result is None:
            continue
        
        summary_data.append({
            'Condition': result['condition_name'],
            'N_Tracks': result['n_tracks'],
            'Mean_Lifetime_s': np.mean(result['lifetimes']),
            'Std_Lifetime_s': np.std(result['lifetimes']),
            'Mean_EndToEnd_um': np.mean(result['end_to_end']),
            'Std_EndToEnd_um': np.std(result['end_to_end']),
            'Mean_PathLength_um': np.mean(result['path_lengths']),
            'Std_PathLength_um': np.std(result['path_lengths']),
            'Mean_FreeExitTime_s': np.mean(result['free_exit_times']) if result['free_exit_times'] else np.nan,
            'Std_FreeExitTime_s': np.std(result['free_exit_times']) if result['free_exit_times'] else np.nan
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(os.path.join(output_dir, 'summary_statistics.csv'), index=False)
    print(f"\nSummary saved to: {os.path.join(output_dir, 'summary_statistics.csv')}")
    print("\n" + "="*60)
    print(summary_df.to_string(index=False))
    print("="*60)

def main():
    """
    Main execution function
    """
    print("\n" + "="*60)
    print("BATCH SMT ANALYSIS - Multiple Conditions")
    print("="*60)
    
    # Get user inputs
    base_dir = input("\nEnter the base directory containing condition folders: ").strip().replace('"', '').replace("'", "")
    output_dir = input("Enter the output directory for results: ").strip().replace('"', '').replace("'", "")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Define conditions and colors
    conditions = {
        'MoeS': 'gold',
        'MoeL': 'orange',
        'MoeSActin': 'red',
        'MoeLActin': 'purple'
    }
    
    # Define radii for FPT analysis
    radii = [1, 2, 3]
    
    # Process each condition
    all_results = []
    for condition_name, color in conditions.items():
        condition_path = os.path.join(base_dir, condition_name)
        
        if not os.path.exists(condition_path):
            print(f"\nWARNING: Directory not found: {condition_path}")
            all_results.append(None)
            continue
        
        result = process_condition(condition_path, condition_name, conditions, radii)
        all_results.append(result)
    
    # Generate all plots
    print("\n" + "="*60)
    print("GENERATING PLOTS")
    print("="*60)
    
    plot_fpt_distributions(all_results, radii, output_dir)
    plot_lifetime_distributions(all_results, output_dir)
    plot_end_to_end_distributions(all_results, output_dir)
    plot_path_length_distributions(all_results, output_dir)
    plot_free_exit_time_distributions(all_results, output_dir)
    
    # Generate summary report
    generate_summary_report(all_results, output_dir)
    
    print("\n" + "="*60)
    print("ANALYSIS COMPLETE!")
    print(f"Results saved to: {output_dir}")
    print("="*60)

if __name__ == "__main__":
    main()
