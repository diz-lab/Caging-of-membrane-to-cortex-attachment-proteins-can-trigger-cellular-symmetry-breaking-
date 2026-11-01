import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.mixture import GaussianMixture
from pathlib import Path
import os
import glob
from datetime import datetime
import warnings
import traceback

# Filter warnings for cleaner output (e.g., from GMM fitting)
warnings.filterwarnings('ignore')

# Set plotting style
sns.set_theme(style="whitegrid")

class UnifiedTrajectoryAnalyzer:
    """
    Complete trajectory analysis from raw spots/tracks to classification.
    
    Workflow:
    1. Load raw spots and tracks files
    2. Compute trajectory features (MBCR, MBCC-CoM, MSD, etc.)
    3. Analyze distributions to suggest thresholds
    4. Classify trajectories using proper two-step logic
    5. Generate comprehensive outputs and visualizations
    """
    
    def __init__(self, input_folder, output_folder, frame_interval=0.02, 
                 pixel_to_micron=0.097, min_points_threshold=5, 
                 min_track_duration=0.2):
        """
        Parameters:
        -----------
        input_folder : str
            Folder containing *-spots.csv and *-tracks.csv files
        output_folder : str
            Folder for all output files
        frame_interval : float
            Time between frames in seconds (default: 0.02 = 20 ms)
        pixel_to_micron : float
            Conversion factor (default: 0.097 µm/pixel)
        min_points_threshold : int
            Minimum localizations required per trajectory for FEATURE CALCULATION (default: 5)
        min_track_duration : float
            Minimum trajectory duration in seconds for CLASSIFICATION (vs 'Short') (default: 0.2 s)
        """
        self.input_folder = Path(input_folder)
        self.output_folder = Path(output_folder)
        self.frame_interval = frame_interval
        self.pixel_to_micron = pixel_to_micron
        self.min_points_threshold = min_points_threshold
        self.min_track_duration = min_track_duration  # in seconds
        
        # Calculate equivalent number of points from time duration
        self.n_points_classification_threshold = int(np.ceil(min_track_duration / frame_interval))
        
        # Classification thresholds (will be refined during analysis)
        self.mbcr_threshold = 0.118  # µm (118 nm) - literature value
        self.mbcc_com_threshold = 0.35 # 35% - literature value
        
        # Create output directory
        self.output_folder.mkdir(parents=True, exist_ok=True)
        
        print("="*80)
        print("UNIFIED TRAJECTORY ANALYZER - RAW DATA TO CLASSIFICATION")
        print("="*80)
        print(f"Input folder: {self.input_folder}")
        print(f"Output folder: {self.output_folder}")
        print(f"Frame interval: {self.frame_interval} s ({self.frame_interval*1000} ms)")
        print(f"Pixel to micron: 1 pixel = {self.pixel_to_micron} µm")
        print(f"Minimum localizations for feature calculation: {self.min_points_threshold}")
        print(f"Minimum track duration for classification: {self.min_track_duration} s")
        print(f"  → Equivalent to ≥ {self.n_points_classification_threshold} points")
        print("="*80)
    
    def find_csv_pairs(self):
        """Find matching spots and tracks files."""
        # Finds files ending with '-spots.csv' or '_spots.csv' and similar for tracks
        spots_files = list(self.input_folder.glob('*-spots.csv')) + list(self.input_folder.glob('*_spots.csv'))
        tracks_files = list(self.input_folder.glob('*-tracks.csv')) + list(self.input_folder.glob('*_tracks.csv'))
        
        pairs = []
        for spots_file in spots_files:
            # Normalize base name (e.g., remove '-spots' and '_spots')
            base_name = spots_file.stem.replace('-spots', '').replace('_spots', '')
            
            for tracks_file in tracks_files:
                if base_name in tracks_file.stem:
                    pairs.append({
                        'base_name': base_name,
                        'spots_file': spots_file,
                        'tracks_file': tracks_file
                    })
                    break
        
        return pairs

    def load_cleaned_csv(self, file_path):
        """Robustly load CSV files, handling extra headers and metadata."""
        # Read the file, skipping rows that might be extra headers
        df = pd.read_csv(file_path)
        # If the first row contains column names (common in some exports), skip it
        if not df.empty and any(df.iloc[0].apply(lambda x: str(x).strip() in df.columns)):
            df = pd.read_csv(file_path, skiprows=1)
        
        # Remove any rows that are identical to the header (extra headers in file)
        for col in df.columns:
            df = df[df[col].astype(str).str.strip() != str(col).strip()]
        
        # Convert numeric columns to float, forcing errors to NaN
        numeric_cols = ['POSITION_X', 'POSITION_Y', 'FRAME', 'TRACK_ID', 
                        'TRACK_DURATION', 'TRACK_MEAN_SPEED', 'TRACK_DISPLACEMENT']
        for col in numeric_cols:
            if col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')
        
        # Drop rows where essential data is NaN after conversion
        if 'TRACK_ID' in df.columns:
            df = df.dropna(subset=['TRACK_ID'])
        
        return df
    
    def calculate_msd_and_diffusion(self, x, y, frames):
        """Calculate MSD and diffusion coefficient with improved fitting."""
        n = len(x)
        # Check if enough points exist for MSD (need at least 3 for a meaningful fit to 1st 4 lags)
        if n < 3:
            return {
                'msd_values': [],
                'time_lags': [],
                'diffusion_coeff': np.nan,
                'msd_at_lag1': np.nan,
                'fit_r_squared': np.nan
            }
        
        # Calculate maximum lag (use up to 1/3 of the track length for robust fitting)
        max_lag = min(n - 1, max(3, int(n / 3)))
        
        msd_values = []
        time_lags = []
        
        # Calculate MSD for each time lag
        for lag in range(1, max_lag + 1):
            squared_displacements = []
            for i in range(n - lag):
                dx = x[i + lag] - x[i]
                dy = y[i + lag] - y[i]
                squared_displacements.append(dx**2 + dy**2)
            
            if len(squared_displacements) > 0:
                msd = np.mean(squared_displacements)
                msd_values.append(msd)
                time_lags.append(lag * self.frame_interval)
        
        if len(time_lags) < 3:
            return {
                'msd_values': msd_values,
                'time_lags': time_lags,
                'diffusion_coeff': np.nan,
                'msd_at_lag1': msd_values[0] if len(msd_values) > 0 else np.nan,
                'fit_r_squared': np.nan
            }
        
        # Fit linear regression to first 4 points (standard practice for local diffusivity)
        fit_points = min(4, len(time_lags))
        t_fit = np.array(time_lags[:fit_points])
        msd_fit = np.array(msd_values[:fit_points])
        
        ss_tot = np.sum((msd_fit - np.mean(msd_fit))**2)
        
        # 1. Fit through origin (MSD(0) = 0)
        if np.sum(t_fit**2) > 0:
            slope_origin = np.sum(t_fit * msd_fit) / np.sum(t_fit**2)
            msd_pred_origin = slope_origin * t_fit
            ss_res_origin = np.sum((msd_fit - msd_pred_origin)**2)
            r_squared_origin = 1 - (ss_res_origin / ss_tot) if ss_tot > 0 else 0
            diffusion_origin = slope_origin / 4.0
        else:
            r_squared_origin = -1 # Assign low score
            diffusion_origin = 0
        
        # 2. Standard linear fit (allows intercept, captures localization error)
        try:
            coeffs = np.polyfit(t_fit, msd_fit, 1)
            slope_intercept = coeffs[0]
            msd_pred_intercept = np.polyval(coeffs, t_fit)
            ss_res_intercept = np.sum((msd_fit - msd_pred_intercept)**2)
            r_squared_intercept = 1 - (ss_res_intercept / ss_tot) if ss_tot > 0 else 0
            diffusion_intercept = slope_intercept / 4.0
        except np.linalg.LinAlgError:
            r_squared_intercept = -1
            diffusion_intercept = 0
        # Use the method with better R²
        if r_squared_origin > r_squared_intercept:
            diffusion_coeff = diffusion_origin
            r_squared = r_squared_origin
        else:
            diffusion_coeff = diffusion_intercept
            r_squared = r_squared_intercept
        
        # Handle negative diffusion coefficients (can happen due to noise/fit choice)
        if diffusion_coeff < 0:
            diffusion_coeff = abs(diffusion_coeff) 
        
        return {
            'msd_values': msd_values,
            'time_lags': time_lags,
            'diffusion_coeff': diffusion_coeff,
            'msd_at_lag1': msd_values[0] if len(msd_values) > 0 else np.nan,
            'fit_r_squared': r_squared
        }
    
    def compute_trajectory_features(self, spots, tracks):
        """
        Compute all trajectory features including MBCR, MBCC-CoM, MSD.
        """
        features = []
        total_tracks = len(tracks)
        
        print(f"  Computing features for {total_tracks} tracks...")
        for idx, track in tracks.iterrows():
            if (idx + 1) % 50 == 0:
                print(f"   Processing track {idx + 1}/{total_tracks}...", end='\r')
            
            track_id = track['TRACK_ID']
            track_spots = spots[spots['TRACK_ID'] == track_id].copy()
            
            # --- CRITICAL FIX 1: CONVERT FRAMES TO SECONDS ---
            raw_duration_frames = float(track['TRACK_DURATION'])
            track_duration_seconds = raw_duration_frames * self.frame_interval
            
            n_points = len(track_spots)
            
            # --- CRITICAL FIX 2: HANDLE "TOO FEW POINTS" TRACKS ---
            if n_points < self.min_points_threshold:
                features.append({
                    'TRACK_ID': track_id,
                    'n_points': n_points,
                    'track_duration': track_duration_seconds,
                    'mean_speed': float(track['TRACK_MEAN_SPEED']) * self.pixel_to_micron,
                    'classification': 'Short', # Pre-classify as Short due to low N
                    # Fill others with NaN
                    'mbcr': np.nan,
                    'mbcc_com_ratio': np.nan,
                    'mbcc_com_distance': np.nan,
                    'total_length': np.nan,
                    'displacement': float(track['TRACK_DISPLACEMENT']) * self.pixel_to_micron,
                    'tortuosity': np.nan,
                    'elongation': np.nan,
                    'msd_at_lag1': np.nan,
                    'diffusion_coeff': np.nan,
                    'msd_fit_r_squared': np.nan,
                    'centroid_x': np.nan,
                    'centroid_y': np.nan,
                    'mbc_center_x': np.nan,
                    'mbc_center_y': np.nan,
                })
                continue
            # ----------------------------------------------------------------
            # 2. PROCESS VALID TRACKS (>= min_points_threshold)
            # ----------------------------------------------------------------
            track_spots = track_spots.sort_values('FRAME')
            
            # Convert pixels to microns
            x = track_spots['POSITION_X'].values * self.pixel_to_micron
            y = track_spots['POSITION_Y'].values * self.pixel_to_micron
            frames = track_spots['FRAME'].values
            
            # Centroid (center of mass)
            centroid_x = np.mean(x)
            centroid_y = np.mean(y)
            
            # Minimum Bounding Circle Radius (MBCR)
            min_max_dist = float('inf')
            best_center = (centroid_x, centroid_y)
            
            candidates = [(centroid_x, centroid_y)] + list(zip(x, y))
            
            for cx, cy in candidates:
                max_dist = np.max(np.sqrt((x - cx)**2 + (y - cy)**2))
                if max_dist < min_max_dist:
                    min_max_dist = max_dist
                    best_center = (cx, cy)
            
            mbcr = min_max_dist
            mbc_center_x, mbc_center_y = best_center
            
            # MBCC-CoM
            mbcc_com_distance = np.sqrt((mbc_center_x - centroid_x)**2 + 
                                       (mbc_center_y - centroid_y)**2)
            mbcc_com_ratio = mbcc_com_distance / mbcr if mbcr > 0 else 0
            
            # Path length and displacement
            displacements = np.sqrt(np.diff(x)**2 + np.diff(y)**2)
            total_length = np.sum(displacements)
            straight_line_distance = np.sqrt((x[-1] - x[0])**2 + (y[-1] - y[0])**2)
            
            # Tortuosity
            if straight_line_distance > 0:
                tortuosity = total_length / straight_line_distance
            else:
                tortuosity = 1.0 if total_length == 0 else float('inf')
            
            # Elongation (Max pairwise distance / MBC Diameter)
            max_pairwise_dist = np.max([np.sqrt((x[i] - x[j])**2 + (y[i] - y[j])**2) 
                                        for i in range(n_points) for j in range(i+1, n_points)])
            elongation = max_pairwise_dist / (2 * mbcr) if mbcr > 0 else 1.0
            # MSD and Diffusion
            msd_results = self.calculate_msd_and_diffusion(x, y, frames)
            
            # CHANGE 3: Check for invalid parameters (NaN or very low values)
            is_invalid = False
            if (np.isnan(mbcr) or mbcr < 0.001 or
                np.isnan(mbcc_com_ratio) or 
                np.isnan(msd_results['diffusion_coeff']) or
                np.isnan(msd_results['msd_at_lag1']) or
                np.isinf(tortuosity)):
                is_invalid = True
            
            features.append({
                'TRACK_ID': track_id,
                'n_points': n_points,
                'mbcr': mbcr,
                'mbcc_com_ratio': mbcc_com_ratio,
                'mbcc_com_distance': mbcc_com_distance,
                'total_length': total_length,
                'displacement': float(track['TRACK_DISPLACEMENT']) * self.pixel_to_micron,
                'tortuosity': tortuosity,
                'elongation': elongation,
                'track_duration': track_duration_seconds, 
                'mean_speed': float(track['TRACK_MEAN_SPEED']) * self.pixel_to_micron, 
                'centroid_x': centroid_x,
                'centroid_y': centroid_y,
                'mbc_center_x': mbc_center_x,
                'mbc_center_y': mbc_center_y,
                'msd_at_lag1': msd_results['msd_at_lag1'],
                'diffusion_coeff': msd_results['diffusion_coeff'],
                'msd_fit_r_squared': msd_results['fit_r_squared'],
                'classification': 'Invalid' if is_invalid else None  # CHANGE 3
            })
        
        print(f"   Processing track {total_tracks}/{total_tracks}... Done!")
        
        return pd.DataFrame(features)
    
    def analyze_thresholds_gmm(self, features_df):
        """
        Analyze distributions and suggest optimal thresholds using GMM.
        """
        print("\n" + "="*80)
        print("THRESHOLD ANALYSIS - FINDING OPTIMAL CUTOFFS")
        print("="*80)
        
        # Filter out tracks too short for feature calculation AND invalid tracks
        analyzable_df = features_df[(features_df['mbcr'].notna()) & 
                                    (features_df['classification'] != 'Invalid')].copy()
        
        # Step 1: Analyze MBCR distribution (all analyzable trajectories)
        mbcr_data = analyzable_df['mbcr'].dropna()
        
        print("\nSTEP 1: MBCR THRESHOLD (Immobile vs Mobile+Hybrid)")
        print("-"*80)
        
        suggested_mbcr = self.mbcr_threshold # Default to literature
        
        if len(mbcr_data) >= 20:
            try:
                X = mbcr_data.values.reshape(-1, 1)
                gmm = GaussianMixture(n_components=2, random_state=42, max_iter=200)
                gmm.fit(X)
                
                means = gmm.means_.flatten()
                sorted_means = np.sort(means)
                
                # Find valley between peaks
                x = np.linspace(sorted_means[0], sorted_means[1], 1000).reshape(-1, 1)
                log_probs = gmm.score_samples(x)
                min_idx = np.argmin(log_probs)
                mbcr_gmm = x[min_idx][0]
                
                print(f"  GMM Peak 1 (Immobile): {sorted_means[0]:.4f} µm")
                print(f"  GMM Peak 2 (Mobile+Hybrid): {sorted_means[1]:.4f} µm")
                print(f"  GMM Valley (suggested threshold): {mbcr_gmm:.4f} µm")
                
                # Use GMM if reasonable (e.g., between 50 nm and 300 nm)
                if 0.05 < mbcr_gmm < 0.3:
                    suggested_mbcr = mbcr_gmm
                    print(f"  → USING GMM threshold: {suggested_mbcr:.4f} µm")
                else:
                    print(f"  → GMM unreasonable, USING literature: {suggested_mbcr:.4f} µm")
            
            except Exception as e:
                print(f"  GMM failed: {e}")
                print(f"  → USING literature value: {suggested_mbcr:.4f} µm")
        else:
            print(f"  Insufficient analyzable data (n={len(mbcr_data)}) for MBCR GMM")
            print(f"  → USING literature value: {suggested_mbcr:.4f} µm")
        
        # Step 2: Analyze MBCC-CoM distribution (only high MBCR trajectories)
        high_mbcr_df = analyzable_df[analyzable_df['mbcr'] >= suggested_mbcr].copy()
        mbcc_com_data = high_mbcr_df['mbcc_com_ratio'].dropna()
        
        print(f"\nSTEP 2: MBCC-CoM THRESHOLD (Mobile vs Hybrid, n={len(mbcc_com_data)} high MBCR trajectories)")
        print("-"*80)
        
        suggested_mbcc = self.mbcc_com_threshold # Default to literature
        
        if len(mbcc_com_data) >= 20:
            try:
                X = mbcc_com_data.values.reshape(-1, 1)
                gmm = GaussianMixture(n_components=2, random_state=42, max_iter=200)
                gmm.fit(X)
                
                means = gmm.means_.flatten()
                sorted_means = np.sort(means)
                
                # Find valley
                x = np.linspace(sorted_means[0], sorted_means[1], 1000).reshape(-1, 1)
                log_probs = gmm.score_samples(x)
                min_idx = np.argmin(log_probs)
                mbcc_gmm = x[min_idx][0]
                
                print(f"  GMM Peak 1 (Mobile): {sorted_means[0]:.4f}")
                print(f"  GMM Peak 2 (Hybrid): {sorted_means[1]:.4f}")
                print(f"  GMM Valley (suggested threshold): {mbcc_gmm:.4f}")
                
                # Use GMM if reasonable (e.g., between 20% and 80%)
                if 0.2 < mbcc_gmm < 0.8:
                    suggested_mbcc = mbcc_gmm
                    print(f"  → USING GMM threshold: {suggested_mbcc:.4f}")
                else:
                    print(f"  → GMM unreasonable, USING literature: {suggested_mbcc:.4f}")
            
            except Exception as e:
                print(f"  GMM failed: {e}")
                print(f"  → USING literature value: {suggested_mbcc:.4f}")
        else:
            print(f"  Insufficient high-MBCR data (n={len(mbcc_com_data)}) for MBCC-CoM GMM")
            print(f"  → USING literature value: {suggested_mbcc:.4f}")
        
        print("="*80)
        
        return suggested_mbcr, suggested_mbcc
    
    def classify_trajectories(self, features_df, mbcr_thresh, mbcc_com_thresh):
        """
        Classify trajectories using proper two-step logic.
        """
        classifications = []
        
        for idx, row in features_df.iterrows():
            
            # Check for pre-classification (Invalid or too few points for feature calculation)
            if row['classification'] == 'Invalid':
                classifications.append('Invalid')
                continue
            elif row['classification'] == 'Short':
                classifications.append('Short')
                continue
            
            track_duration = row['track_duration'] # Already in seconds
            mbcr = row['mbcr'] 
            mbcc_com = row['mbcc_com_ratio']
            
            # 1. Check by duration 
            if track_duration < self.min_track_duration:
                classification = 'Short'
            # 2. Step 1: MBCR check
            elif mbcr < mbcr_thresh:
                classification = 'Immobile'
            # 3. Step 2: MBCC-CoM check (for high MBCR only)
            else:
                if mbcc_com < mbcc_com_thresh:
                    classification = 'Mobile'
                else:
                    classification = 'Hybrid'
            
            classifications.append(classification)
        
        features_df['classification'] = classifications
        return features_df
    
    def create_visualizations(self, features_df, base_name, mbcr_thresh, mbcc_thresh, 
                            suggested_mbcr, suggested_mbcc):
        """
        Create comprehensive visualization suite.
        """
        colors = {'Mobile': '#C41E3A', 'Immobile': '#1E40AF', 'Hybrid': '#D4AF37', 
                 'Short': '#808080', 'Invalid': '#FF0000'}
        
        fig = plt.figure(figsize=(20, 12))
        gs = fig.add_gridspec(3, 4, hspace=0.3, wspace=0.3)
        
        # Filter out Invalid tracks for plotting
        plot_df = features_df[features_df['classification'] != 'Invalid'].copy()
        
        # 1. Classification Distribution (Pie chart: Mobile, Hybrid, Immobile only)
        ax = fig.add_subplot(gs[0, 0])
        analyzed_categories = ['Mobile', 'Hybrid', 'Immobile']
        class_counts = features_df[features_df['classification'].isin(analyzed_categories)]['classification'].value_counts()
        if not class_counts.empty:
            ax.pie(class_counts.values, labels=class_counts.index, autopct='%1.1f%%',
                   colors=[colors.get(c, 'gray') for c in class_counts.index], startangle=90)
        ax.set_title('Analyzed Tracks Distribution\n(Mobile/Hybrid/Immobile)', fontsize=14, fontweight='bold')
        
        # 2. MBCR Distribution
        ax = fig.add_subplot(gs[0, 1])
        for classification in ['Immobile', 'Mobile', 'Hybrid']:
            data = plot_df[plot_df['classification'] == classification]['mbcr'].dropna()
            if len(data) > 0:
                ax.hist(data, bins=30, alpha=0.6, label=classification, color=colors[classification])
        ax.axvline(suggested_mbcr, color='gray', linestyle=':', linewidth=2, label=f'Suggested ({suggested_mbcr:.3f})')
        ax.axvline(mbcr_thresh, color='black', linestyle='-', linewidth=2, label=f'Applied ({mbcr_thresh:.3f})')
        ax.set_xlabel('MBCR (µm)', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title('MBCR Distribution', fontsize=14, fontweight='bold')
        ax.legend()
        
        # 3. MBCC-CoM Distribution
        ax = fig.add_subplot(gs[0, 2])
        for classification in ['Mobile', 'Hybrid']:
            data = plot_df[plot_df['classification'] == classification]['mbcc_com_ratio'].dropna()
            if len(data) > 0:
                ax.hist(data, bins=30, alpha=0.6, label=classification, color=colors[classification])
        ax.axvline(suggested_mbcc, color='gray', linestyle=':', linewidth=2, label=f'Suggested ({suggested_mbcc:.2f})')
        ax.axvline(mbcc_thresh, color='black', linestyle='-', linewidth=2, label=f'Applied ({mbcc_thresh:.2f})')
        ax.set_xlabel('MBCC-CoM Ratio', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title('MBCC-CoM Distribution', fontsize=14, fontweight='bold')
        ax.legend()
        
        # 4. Decision Space
        ax = fig.add_subplot(gs[0, 3])
        for classification in ['Immobile', 'Mobile', 'Hybrid']:
            data = plot_df[plot_df['classification'] == classification].dropna(subset=['mbcr', 'mbcc_com_ratio'])
            if len(data) > 0:
                ax.scatter(data['mbcr'], data['mbcc_com_ratio'], alpha=0.5, s=20, label=classification, color=colors[classification])
        ax.axvline(suggested_mbcr, color='gray', linestyle=':', linewidth=1.5, alpha=0.5)
        ax.axvline(mbcr_thresh, color='black', linestyle='-', linewidth=1.5, alpha=0.7)
        ax.axhline(suggested_mbcc, color='gray', linestyle=':', linewidth=1.5, alpha=0.5)
        ax.axhline(mbcc_thresh, color='black', linestyle='-', linewidth=1.5, alpha=0.7)
        ax.set_xlabel('MBCR (µm)', fontsize=12)
        ax.set_ylabel('MBCC-CoM Ratio', fontsize=12)
        ax.set_title('Decision Space', fontsize=14, fontweight='bold')
        ax.legend()
        
        # 5. Diffusion Coefficient
        ax = fig.add_subplot(gs[1, 0])
        for classification in ['Immobile', 'Mobile', 'Hybrid']:
            data = plot_df[plot_df['classification'] == classification]['diffusion_coeff'].dropna()
            if len(data) > 0:
                ax.hist(data, bins=30, alpha=0.6, label=classification, color=colors[classification])
        ax.set_xlabel('Diffusion Coeff (µm²/s)', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title('Diffusion Distribution', fontsize=14, fontweight='bold')
        ax.legend()
        
        # 6. MSD at lag 1
        ax = fig.add_subplot(gs[1, 1])
        for classification in ['Immobile', 'Mobile', 'Hybrid']:
            data = plot_df[plot_df['classification'] == classification]['msd_at_lag1'].dropna()
            if len(data) > 0:
                ax.hist(data, bins=30, alpha=0.6, label=classification, color=colors[classification])
        ax.set_xlabel('MSD at lag=1 (µm²)', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title('MSD Distribution', fontsize=14, fontweight='bold')
        ax.legend()
        
        # 7. Tortuosity
        ax = fig.add_subplot(gs[1, 2])
        for classification in ['Immobile', 'Mobile', 'Hybrid']:
            data = plot_df[plot_df['classification'] == classification]['tortuosity'].dropna()
            data = data[data < 10] # Filter outliers for visualization
            if len(data) > 0:
                ax.hist(data, bins=30, alpha=0.6, label=classification, color=colors[classification])
        ax.set_xlabel('Tortuosity', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title('Tortuosity Distribution', fontsize=14, fontweight='bold')
        ax.legend()
        
        # 8. Track Duration
        ax = fig.add_subplot(gs[1, 3])
        for classification in ['Immobile', 'Mobile', 'Hybrid', 'Short']:
            data = plot_df[plot_df['classification'] == classification]['track_duration'].dropna()
            if len(data) > 0:
                ax.hist(data, bins=30, alpha=0.6, label=classification, color=colors[classification])
        ax.axvline(self.min_track_duration, color='black', linestyle='--', linewidth=2, label=f'Min dur ({self.min_track_duration}s)')
        ax.set_xlabel('Track Duration (s)', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title('Track Duration Distribution', fontsize=14, fontweight='bold')
        ax.legend()
        
        plt.suptitle(f'Trajectory Analysis: {base_name}', fontsize=18, fontweight='bold', y=0.98)
        output_file = self.output_folder / f'{base_name}_analysis.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        return output_file

    def create_text_report(self, features_df, base_name, mbcr_thresh, mbcc_thresh):
        """Create detailed text report."""
        text_file = self.output_folder / f'{base_name}_report.txt'
        n_analyzed = len(features_df[~features_df['classification'].isin(['Short', 'Invalid'])])
        
        with open(text_file, 'w') as f:
            f.write("="*80 + "\n")
            f.write(f"UNIFIED TRAJECTORY ANALYSIS REPORT: {base_name}\n")
            f.write("="*80 + "\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Frame interval: {self.frame_interval} s ({self.frame_interval*1000} ms)\n")
            f.write(f"Pixel to micron: 1 pixel = {self.pixel_to_micron} µm\n")
            f.write(f"Min points for features: {self.min_points_threshold}\n")
            f.write(f"Min duration for classification: {self.min_track_duration} s\n")
            f.write("="*80 + "\n\n")
            
            f.write(f"MBCR threshold: {mbcr_thresh:.4f} µm\n")
            f.write(f"MBCC-CoM threshold: {mbcc_thresh:.4f}\n\n")
            
            f.write("CLASSIFICATION RESULTS\n")
            f.write("-"*80 + "\n")
            
            for classification in ['Immobile', 'Mobile', 'Hybrid']:
                class_data = features_df[features_df['classification'] == classification]
                count = len(class_data)
                if n_analyzed > 0:
                    pct = (count / n_analyzed) * 100
                else:
                    pct = 0
                f.write(f"{classification}: {count} ({pct:.1f}% of analyzed)\n")
                
                if count > 0:
                    f.write(f"  Mean MBCR: {class_data['mbcr'].mean():.4f} µm\n")
                    f.write(f"  Mean MBCC-CoM: {class_data['mbcc_com_ratio'].mean():.4f}\n")
                    f.write(f"  Mean Diffusion Coeff: {class_data['diffusion_coeff'].mean():.4f} µm²/s\n")
                    f.write(f"  Mean MSD (lag=1): {class_data['msd_at_lag1'].mean():.4f} µm²\n")
                    f.write(f"  Mean Tortuosity: {class_data['tortuosity'].mean():.4f}\n")
                    f.write(f"  Mean Speed: {class_data['mean_speed'].mean():.4f} µm/s\n\n")
            
            f.write("="*80 + "\n")
        
        return text_file

    def process_file_pair(self, pair, suggested_mbcr=None, suggested_mbcc=None):
        """
        Process a single spots/tracks file pair.
        """
        base_name = pair['base_name']
        spots_file = pair['spots_file']
        tracks_file = pair['tracks_file']
        
        print("\n" + "="*80)
        print(f"PROCESSING: {base_name}")
        print("="*80)
        print(f"Spots file: {spots_file.name}")
        print(f"Tracks file: {tracks_file.name}")
        
        try:
            # Load data with robust header handling
            print("\nLoading data...")
            spots = self.load_cleaned_csv(spots_file)
            tracks = self.load_cleaned_csv(tracks_file)
            print(f"  Loaded {len(spots)} spots and {len(tracks)} tracks")
            
            # Compute features
            print("\nComputing trajectory features...")
            features_df = self.compute_trajectory_features(spots, tracks)
            
            # Analyze thresholds (only if not provided from batch)
            if suggested_mbcr is None or suggested_mbcc is None:
                suggested_mbcr, suggested_mbcc = self.analyze_thresholds_gmm(features_df)
            
            # Use suggested thresholds as applied thresholds
            mbcr_thresh = suggested_mbcr
            mbcc_thresh = suggested_mbcc
            
            # Classify trajectories
            print("\nClassifying trajectories...")
            features_df = self.classify_trajectories(features_df, mbcr_thresh, mbcc_thresh)
            
            # Create visualizations
            print("\nCreating visualizations...")
            viz_file = self.create_visualizations(features_df, base_name, mbcr_thresh, mbcc_thresh,
                                                 suggested_mbcr, suggested_mbcc)
            print(f"  Saved: {viz_file.name}")
            
            # Save results
            print("\nSaving results...")
            csv_comp, csv_summ, excel = self.save_results(features_df, base_name)
            print(f"  Saved: {csv_comp.name}")
            print(f"  Saved: {csv_summ.name}")
            print(f"  Saved: {excel.name}")
            
            # Create text report
            print("\nCreating text report...")
            report_file = self.create_text_report(features_df, base_name, mbcr_thresh, mbcc_thresh)
            print(f"  Saved: {report_file.name}")
            
            print(f"\n✓ Successfully processed {base_name}")
            
            return features_df, suggested_mbcr, suggested_mbcc, mbcr_thresh, mbcc_thresh
            
        except Exception as e:
            print(f"\n✗ Error processing {base_name}: {str(e)}")
            traceback.print_exc()
            return None, None, None, None, None

    def create_master_visualization(self, all_features, mbcr_thresh, mbcc_thresh,
                                   suggested_mbcr, suggested_mbcc):
        """
        Create master visualization combining all FOVs.
        """
        colors = {'Mobile': '#C41E3A', 'Immobile': '#1E40AF', 'Hybrid': '#D4AF37', 
                 'Short': '#808080', 'Invalid': '#FF0000'}
        
        fig = plt.figure(figsize=(24, 14))
        gs = fig.add_gridspec(3, 5, hspace=0.3, wspace=0.3)
        
        # Filter data
        plot_df = all_features[all_features['classification'] != 'Invalid'].copy()
        analyzed_df = all_features[all_features['classification'].isin(['Mobile', 'Hybrid', 'Immobile'])]
        
        # 1. Overall Classification Distribution (Pie chart: Mobile, Hybrid, Immobile only)
        ax = fig.add_subplot(gs[0, 0])
        class_counts_analyzed = analyzed_df['classification'].value_counts()
        if not class_counts_analyzed.empty:
            ax.pie(class_counts_analyzed.values, labels=class_counts_analyzed.index, autopct='%1.1f%%',
                   colors=[colors.get(c, 'gray') for c in class_counts_analyzed.index], startangle=90)
        ax.set_title('Analyzed Tracks Only\n(Mobile/Hybrid/Immobile)', 
                    fontsize=13, fontweight='bold')
        
        # 2. MBCR Distribution
        ax = fig.add_subplot(gs[0, 2])
        for classification in ['Immobile', 'Mobile', 'Hybrid']:
            data = plot_df[plot_df['classification'] == classification]['mbcr'].dropna()
            if len(data) > 0:
                ax.hist(data, bins=50, alpha=0.6, label=f"{classification} (n={len(data)})", 
                       color=colors[classification])
        ax.axvline(suggested_mbcr, color='gray', linestyle=':', linewidth=2.5, 
                  label=f'Suggested ({suggested_mbcr:.3f})')
        ax.axvline(mbcr_thresh, color='black', linestyle='-', linewidth=2.5, 
                  label=f'Applied ({mbcr_thresh:.3f})')
        ax.set_xlabel('MBCR (µm)', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title('Master MBCR Distribution', fontsize=14, fontweight='bold')
        ax.legend()
        
        # 3. MBCC-CoM Distribution
        ax = fig.add_subplot(gs[0, 3])
        for classification in ['Mobile', 'Hybrid']:
            data = plot_df[plot_df['classification'] == classification]['mbcc_com_ratio'].dropna()
            if len(data) > 0:
                ax.hist(data, bins=50, alpha=0.6, label=f"{classification} (n={len(data)})",
                       color=colors[classification])
        ax.axvline(suggested_mbcc, color='gray', linestyle=':', linewidth=2.5,
                  label=f'Suggested ({suggested_mbcc:.2f})')
        ax.axvline(mbcc_thresh, color='black', linestyle='-', linewidth=2.5,
                  label=f'Applied ({mbcc_thresh:.2f})')
        ax.set_xlabel('MBCC-CoM Ratio', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title('Master MBCC-CoM Distribution', fontsize=14, fontweight='bold')
        ax.legend()
        
        # 4. Track Duration
        ax = fig.add_subplot(gs[0, 4])
        for classification in ['Immobile', 'Mobile', 'Hybrid']:
            data = plot_df[plot_df['classification'] == classification]['track_duration'].dropna()
            if len(data) > 0:
                ax.hist(data, bins=50, alpha=0.6, label=classification, 
                       color=colors[classification])
        ax.axvline(self.min_track_duration, color='black', linestyle='--', linewidth=2)
        ax.set_xlabel('Track Duration (s)', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title('Master Track Duration', fontsize=14, fontweight='bold')
        ax.legend()
        
        # 5. Decision Space
        ax = fig.add_subplot(gs[1, 0:2])
        for classification in ['Immobile', 'Mobile', 'Hybrid']:
            data = plot_df[plot_df['classification'] == classification].dropna(subset=['mbcr', 'mbcc_com_ratio'])
            if len(data) > 0:
                ax.scatter(data['mbcr'], data['mbcc_com_ratio'], alpha=0.5, s=20,
                          label=f"{classification} (n={len(data)})", 
                          color=colors[classification])
        ax.axvline(suggested_mbcr, color='gray', linestyle=':', linewidth=2, alpha=0.5)
        ax.axvline(mbcr_thresh, color='black', linestyle='-', linewidth=2, alpha=0.7)
        ax.axhline(suggested_mbcc, color='gray', linestyle=':', linewidth=2, alpha=0.5)
        ax.axhline(mbcc_thresh, color='black', linestyle='-', linewidth=2, alpha=0.7)
        ax.set_xlabel('MBCR (µm)', fontsize=12)
        ax.set_ylabel('MBCC-CoM Ratio', fontsize=12)
        ax.set_title('Master Decision Space', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(alpha=0.3)
        
        # 6. Diffusion Coefficient
        ax = fig.add_subplot(gs[1, 2])
        for classification in ['Immobile', 'Mobile', 'Hybrid']:
            data = plot_df[plot_df['classification'] == classification]['diffusion_coeff'].dropna()
            if len(data) > 0:
                ax.hist(data, bins=50, alpha=0.6, label=classification,
                       color=colors[classification])
        ax.set_xlabel('Diffusion Coeff (µm²/s)', fontsize=12)
        ax.set_ylabel('Count', fontsize=12)
        ax.set_title('Master Diffusion Distribution', fontsize=14, fontweight='bold')
        ax.legend()
        
        # 7. MSD Distribution
        ax = fig.add_subplot(gs[1, 3])
        for classification in ['Immobile', 'Mobile', 'Hybrid']:
            data = all_features[all_features['classification'] == classification]['msd_at_lag1']
            data = data[data.notna()]
            if len(data) > 0:
                ax.hist(data, bins=30, alpha=0.6, label=classification, color=colors[classification])
        ax.set_xlabel('MSD at lag=1 (µm²)', fontsize=11)
        ax.set_ylabel('Count', fontsize=11)
        ax.set_title('MSD Distribution', fontsize=14, fontweight='bold')
        ax.legend()
            
        # 8. Tortuosity Distribution
        ax = fig.add_subplot(gs[2, 0])
        for classification in ['Immobile', 'Mobile', 'Hybrid']:
            data = all_features[all_features['classification'] == classification]['tortuosity']
            data = data[(data.notna()) & (data < 10)]
            if len(data) > 0:
                ax.hist(data, bins=30, alpha=0.6, label=classification, color=colors[classification])
        ax.set_xlabel('Tortuosity', fontsize=11)
        ax.set_ylabel('Count', fontsize=11)
        ax.set_title('Tortuosity Distribution', fontsize=14, fontweight='bold')
        ax.legend()
            
        # 9. Speed Distribution
        ax = fig.add_subplot(gs[2, 1])
        for classification in ['Immobile', 'Mobile', 'Hybrid']:
            data = all_features[all_features['classification'] == classification]['mean_speed']
            if len(data) > 0:
                ax.hist(data, bins=30, alpha=0.6, label=classification, color=colors[classification])
        ax.set_xlabel('Mean Speed (µm/s)', fontsize=11)
        ax.set_ylabel('Count', fontsize=11)
        ax.set_title('Speed Distribution', fontsize=14, fontweight='bold')
        ax.legend()
            
        # 10. Number of Points Distribution
        ax = fig.add_subplot(gs[2, 2])
        for classification in ['Immobile', 'Mobile', 'Hybrid']:
            data = all_features[all_features['classification'] == classification]['n_points']
            if len(data) > 0:
                ax.hist(data, bins=30, alpha=0.6, label=classification, color=colors[classification], log=True)
        ax.axvline(self.min_points_threshold, color='black', linestyle='--', linewidth=2, label=f'Min points ({self.min_points_threshold})')
        ax.set_xlabel('Number of Localizations', fontsize=11)
        ax.set_ylabel('Count (Log Scale)', fontsize=11)
        ax.set_title('Localization Count Distribution', fontsize=14, fontweight='bold')
        ax.legend()
            
        # 11. Statistics Table
        ax = fig.add_subplot(gs[2, 3])
        ax.axis('off')
            
        # Calculate statistics
        total_all = len(all_features)
        n_short = len(all_features[all_features['classification'] == 'Short'])
        n_invalid = len(all_features[all_features['classification'] == 'Invalid'])
        analyzed_total = len(analyzed_df)
        
        stats_text = f"Classification Summary\n{'-'*40}\n"
        stats_text += f"Total Trajectories: {total_all}\n"
        stats_text += f"Short (< {self.min_track_duration}s): {n_short}\n"
        stats_text += f"Invalid (NaN/Low Signal): {n_invalid}\n"
        stats_text += f"Analyzed: {analyzed_total}\n\n"
        stats_text += f"Thresholds Used:\n"
        stats_text += f" MBCR: {mbcr_thresh:.4f} µm\n"
        stats_text += f" MBCC-CoM: {mbcc_thresh:.4f}\n"
        stats_text += f" Min duration: {self.min_track_duration}s\n\n"
        stats_text += f"Distribution (of analyzed):\n"
        
        for classification in ['Immobile', 'Mobile', 'Hybrid']:
            count = len(all_features[all_features['classification'] == classification])
            if analyzed_total > 0:
                pct = (count / analyzed_total) * 100
                stats_text += f"{classification}: {count} ({pct:.1f}%)\n"
            else:
                stats_text += f"{classification}: {count}\n"
                
        ax.text(0.1, 0.9, stats_text, fontsize=10, verticalalignment='top', family='monospace', 
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))
        
        plt.suptitle(f'Comprehensive Trajectory Analysis', fontsize=16, fontweight='bold', y=0.995)
        output_file = self.output_folder / 'master_analysis.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        return output_file

    def save_results(self, features_df, base_name):
        """Save all results to CSV and Excel files."""
        # Comprehensive CSV
        comprehensive_csv = self.output_folder / f'{base_name}_comprehensive_results.csv'
        features_df.to_csv(comprehensive_csv, index=False)
        
        # Summary CSV
        summary_df = features_df[['TRACK_ID', 'classification', 'n_points', 'mbcr', 'mbcc_com_ratio', 'diffusion_coeff', 'msd_at_lag1', 'tortuosity', 'mean_speed', 'track_duration']]
        summary_csv = self.output_folder / f'{base_name}_summary.csv'
        summary_df.to_csv(summary_csv, index=False)
        
        # Excel with multiple sheets
        excel_file = self.output_folder / f'{base_name}_results.xlsx'
        with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
            summary_df.to_excel(writer, sheet_name='Summary', index=False)
            
            # Statistics sheet
            stats_data = []
            for classification in ['Immobile', 'Mobile', 'Hybrid', 'Short', 'Invalid']:
                class_data = features_df[features_df['classification'] == classification]
                if len(class_data) > 0:
                    stats_data.append({
                        'Classification': classification,
                        'Count': len(class_data),
                        'Percentage': f"{len(class_data)/len(features_df)*100:.1f}%",
                        'Mean_MSD': class_data['msd_at_lag1'].mean(),
                        'Mean_Diffusion_Coeff': class_data['diffusion_coeff'].mean(),
                        'Mean_MBCR': class_data['mbcr'].mean(),
                        'Mean_MBCC_CoM': class_data['mbcc_com_ratio'].mean(),
                        'Mean_Tortuosity': class_data['tortuosity'].mean(),
                        'Mean_Speed': class_data['mean_speed'].mean()
                    })
            stats_df = pd.DataFrame(stats_data)
            stats_df.to_excel(writer, sheet_name='Statistics', index=False)
            
            # All features
            features_df.to_excel(writer, sheet_name='All_Features', index=False)
            
        return comprehensive_csv, summary_csv, excel_file

    def generate_master_report(self, all_features_list):
        """Generate a master file and plot of distribution of all tracks across all FOVs."""
        if not all_features_list:
            return
            
        master_df = pd.concat(all_features_list, ignore_index=True)
        master_file = self.output_folder / 'master_comprehensive_results.csv'
        master_df.to_csv(master_file, index=False)
        
        # Master Distribution Plot
        colors = {'Mobile': '#C41E3A', 'Immobile': '#1E40AF', 'Hybrid': '#D4AF37', 
                  'Short': '#808080', 'Invalid': '#000000'}
        
        plt.figure(figsize=(12, 8))
        
        # Prepare data for plotting
        categories = ['All', 'Invalid', 'Analyzed', 'Mobile', 'Hybrid', 'Immobile', 'Short']
        counts = []
        total = len(master_df)
        
        for cat in categories:
            if cat == 'All':
                counts.append(total)
            elif cat == 'Analyzed':
                counts.append(len(master_df[~master_df['classification'].isin(['Short', 'Invalid'])]))
            else:
                counts.append(len(master_df[master_df['classification'] == cat]))
        
        bars = plt.bar(categories, counts, color=['#333333', '#000000', '#555555', '#C41E3A', '#D4AF37', '#1E40AF', '#808080'])
        
        # Add labels
        for bar in bars:
            height = bar.get_height()
            plt.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                     f'{int(height)}\n({height/total*100:.1f}%)',
                     ha='center', va='bottom', fontsize=10)
        
        plt.title('Master Track Distribution (All FOVs)', fontsize=16, fontweight='bold')
        plt.ylabel('Count', fontsize=12)
        plt.grid(axis='y', alpha=0.3)
        
        master_plot = self.output_folder / 'master_distribution_plot.png'
        plt.savefig(master_plot, dpi=300, bbox_inches='tight')
        plt.close()
        
        # Master Text Summary
        master_report = self.output_folder / 'master_report.txt'
        with open(master_report, 'w') as f:
            f.write("="*80 + "\n")
            f.write("MASTER TRAJECTORY ANALYSIS REPORT (ALL FOVs)\n")
            f.write("="*80 + "\n")
            f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Total FOVs processed: {len(all_features_list)}\n")
            f.write(f"Total Trajectories: {total}\n")
            f.write("="*80 + "\n\n")
            
            f.write("Global Classification Summary:\n")
            for cat, count in zip(categories, counts):
                f.write(f"  {cat}: {count} ({count/total*100:.1f}%)\n")
                
        print(f"\nMaster report generated: {master_report}")

    def run_analysis(self):
        """Main execution pipeline with batch processing."""
        pairs = self.find_csv_pairs()
        if not pairs:
            print("No matching spots/tracks files found!")
            return
            
        print(f"Found {len(pairs)} FOV pairs to process.")
        all_features_list = []
        
        for pair in pairs:
            base_name = pair['base_name']
            print(f"\nProcessing: {base_name}")
            
            try:
                # Load data with robust header handling
                spots = self.load_cleaned_csv(pair['spots_file'])
                tracks = self.load_cleaned_csv(pair['tracks_file'])
                
                # Compute features
                features_df = self.compute_trajectory_features(spots, tracks)
                
                # Analyze thresholds
                suggested_mbcr, suggested_mbcc = self.analyze_thresholds_gmm(features_df)
                
                # Use suggested thresholds as applied thresholds (can be manually overridden here)
                mbcr_thresh = suggested_mbcr
                mbcc_thresh = suggested_mbcc
                
                # Classify
                features_df = self.classify_trajectories(features_df, mbcr_thresh, mbcc_thresh)
                
                # Save and visualize
                self.create_visualizations(features_df, base_name, mbcr_thresh, mbcc_thresh, suggested_mbcr, suggested_mbcc)
                self.save_results(features_df, base_name)
                self.create_text_report(features_df, base_name, mbcr_thresh, mbcc_thresh)
                
                # Add to master list
                features_df['FOV'] = base_name
                all_features_list.append(features_df)
                
                print(f"Finished processing {base_name}")
            except Exception as e:
                print(f"Error processing {base_name}: {e}")
                traceback.print_exc()
        
        # Generate master report if multiple FOVs were processed
        if len(all_features_list) > 0:
            # Combine all features for master visualization
            all_features_combined = pd.concat(all_features_list, ignore_index=True)
            # Use the last pair's thresholds for the master plot (or average them)
            self.create_master_visualization(all_features_combined, mbcr_thresh, mbcc_thresh, suggested_mbcr, suggested_mbcc)
            self.generate_master_report(all_features_list)

def get_user_input(prompt, default, type_func=str):
    """Helper function to get user input with a default value."""
    user_val = input(f"{prompt} [{default}]: ").strip()
    if not user_val:
        return default
    try:
        return type_func(user_val)
    except ValueError:
        print(f"Invalid input. Using default: {default}")
        return default

if __name__ == "__main__":
    print("="*80)
    print("UNIFIED TRAJECTORY ANALYZER - CONFIGURATION")
    print("="*80)
    
    input_folder = get_user_input("Input folder", "data/input")
    output_folder = get_user_input("Output folder", "data/output")
    frame_interval = get_user_input("Frame interval (s)", 0.02, float)
    pixel_to_micron = get_user_input("Pixel to micron conversion", 0.097, float)
    min_points = get_user_input("Minimum localizations for features", 5, int)
    min_duration = get_user_input("Minimum track duration for classification (s)", 0.2, float)
    
    analyzer = UnifiedTrajectoryAnalyzer(
        input_folder=input_folder,
        output_folder=output_folder,
        frame_interval=frame_interval,
        pixel_to_micron=pixel_to_micron,
        min_points_threshold=min_points,
        min_track_duration=min_duration
    )
    analyzer.run_analysis()
