#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Grid-based CMD fitting using KDE and Jensen-Shannon divergence
Created on Mon Jul  7 19:29:13 2025
@author: george
"""

import numpy as np
import pandas as pd
import os
import glob
import re
import sys
from KDEpy import FFTKDE
from scipy.spatial.distance import jensenshannon

class GridCMDAnalyzer:
    """
    Streamlined class for comparing observed CMD data with a grid of simulated CMDs
    """
    #    gaia_file = 
    def __init__(self, 
                 run_number,
                 target_age,
                 grid_path=None,
                 #obs_path="/dartfs-hpc/rc/home/8/f007058/cleaned_ngc6791_no_blue_stragglers.csv",
                 obs_path="/dartfs-hpc/rc/home/8/f007058/NGC6791_referee_isochrones_sCMDs/chaboyer_DM13.280_EBV0.200_age8.3Gyr_001.csv",
                 output_path="/dartfs-hpc/rc/home/8/f007058/OUTPUT_KDE",
                 color_range=(1.0, 1.5),
                 mag_range=(15, 19),
                 kde_params={'bw': 0.15, 'grid_size': 1000, 'epsilon': 1e-6}):
        
        self.run_number = run_number
        self.target_age = target_age
        
        # Set grid path based on run number if not provided
        if grid_path is None:
            #self.grid_path = f"/dartfs-hpc/rc/lab/C/ChaboyerB/sCMD_July_11/{run_number}"#/dartfs-hpc/scratch/f007058/sCMD_W
            self.grid_path = f"/dartfs-hpc/scratch/f007058/sCMD_W/{run_number}"
        else:
            self.grid_path = grid_path
            
        self.obs_path = obs_path
        self.output_path = output_path
        self.color_min, self.color_max = color_range
        self.mag_min, self.mag_max = mag_range
        
        # KDE parameters
        self.bw = kde_params['bw']
        self.grid_size = kde_params['grid_size']
        self.epsilon = kde_params['epsilon']
        
        # Create output directory if it doesn't exist
        os.makedirs(self.output_path, exist_ok=True)
        
        # Load and prepare data
        self.obs_data = self._load_observed_data()
        self.grid_data = self._load_grid_data()

        # Precompute eval grid + observed PDF once (reused for every model)
        color_grid = np.linspace(self.color_min - self.epsilon,
                                 self.color_max + self.epsilon,
                                 self.grid_size)
        mag_grid = np.linspace(self.mag_min - self.epsilon,
                               self.mag_max + self.epsilon,
                               self.grid_size)
        color_mesh, mag_mesh = np.meshgrid(color_grid, mag_grid)
        eval_points = np.column_stack([color_mesh.ravel(), mag_mesh.ravel()])
        sort_indices = np.lexsort((eval_points[:, 1], eval_points[:, 0]))
        self._eval_points = eval_points[sort_indices]

        obs_clean = self.obs_data[np.isfinite(self.obs_data).all(axis=1)]
        if len(obs_clean) >= 10:
            kde_obs = FFTKDE(kernel='gaussian', bw=self.bw, norm=1).fit(obs_clean)
            pdf_obs = kde_obs.evaluate(self._eval_points)
            pdf_obs = np.maximum(pdf_obs, 1e-10)
            self._pdf_obs = pdf_obs / np.sum(pdf_obs)
        else:
            self._pdf_obs = None
        
        print(f"✅ Initialized analyzer for run {run_number}, age {target_age:.1f} Gyr")
        print(f"✅ Grid path: {self.grid_path}")
        print(f"✅ Found {len(self.grid_data)} grid models for target age")
        print(f"✅ Observational data: {len(self.obs_data)} stars")
    
    def _load_observed_data(self):
        """Load and filter observational data"""
        print(f"📖 Loading observed data from {self.obs_path}")
        
        data = pd.read_csv(self.obs_path)
        print(f"   Raw data: {len(data)} stars")
        
        # Handle BP-RP color column
        if "BP_RP" not in data.columns:
            if "BPmag" in data.columns and "RPmag" in data.columns:
                data["BP_RP"] = data["BPmag"] - data["RPmag"]
            else:
                raise ValueError("Cannot find BP-RP color data")
        
        # Filter to analysis region
        filtered = data[
            (data["BP_RP"] >= self.color_min) & 
            (data["BP_RP"] <= self.color_max) & 
            (data["Gmag"] >= self.mag_min) & 
            (data["Gmag"] <= self.mag_max)
        ]
        
        print(f"   Filtered data: {len(filtered)} stars")
        
        # Return as numpy array [color, magnitude]
        return filtered[['BP_RP', 'Gmag']].to_numpy()
    
    def _parse_grid_filename(self, filename):
        """Extract parameters from grid filename with updated pattern"""
        # Pattern: NGC6791_age8.0Gyr_DM13.35_EBV0.15chaboyer_DM13.280_EBV0.200_age8.3Gyr_001.csv
        pattern = r'NGC6791_age(\d+\.?\d*)Gyr_DM(\d+\.?\d*)_EBV(\d+\.?\d*)\.csv'
        match = re.search(pattern, filename)
        
        if match:
            return {
                'age': float(match.group(1)),
                'dm': float(match.group(2)),
                'ebv': float(match.group(3))
            }
        else:
            raise ValueError(f"Could not parse filename: {filename}")
    
    def _load_grid_data(self):
        """Load all grid files for the specified age"""
        print(f"📖 Loading grid data from {self.grid_path} for age {self.target_age:.1f} Gyr")
        
        if not os.path.exists(self.grid_path):
            raise FileNotFoundError(f"Grid directory not found: {self.grid_path}")
        
        # Look for files with the specific age pattern
        age_pattern = f"NGC6791_age{self.target_age:.1f}Gyr_*.csv"
        csv_files = glob.glob(f"{self.grid_path}/{age_pattern}")
        
        if not csv_files:
            raise FileNotFoundError(f"No files found matching pattern: {age_pattern} in {self.grid_path}")
        
        print(f"   Found {len(csv_files)} CSV files for age {self.target_age:.1f} Gyr")
        
        grid_data = {}
        loaded_count = 0
        
        for file_path in csv_files:
            filename = os.path.basename(file_path)
            
            try:
                # Parse parameters
                params = self._parse_grid_filename(filename)
                
                # Double-check age matches (should be exact since we filtered by pattern)
                if abs(params['age'] - self.target_age) > 0.01:
                    print(f"   Warning: Age mismatch in {filename}, expected {self.target_age:.1f}, got {params['age']:.1f}")
                    continue
                
                # Load data
                data = pd.read_csv(file_path)
                
                # Store with tuple key for easy access
                key = (params['age'], params['dm'], params['ebv'])
                grid_data[key] = {
                    'data': data,
                    'params': params,
                    'filename': filename,
                    'filepath': file_path
                }
                loaded_count += 1
                
            except Exception as e:
                print(f"   Warning: Skipped {filename} - {e}")
                continue
        
        print(f"   Successfully loaded: {loaded_count} files")
        return grid_data
    
    def _compute_kde_distance(self, sim_color, sim_mag):
        """Compute Jensen-Shannon distance between observed and simulated data using KDE"""

        if self._pdf_obs is None:
            return np.inf

        eval_points = self._eval_points
        pdf_obs = self._pdf_obs

        sim_data = np.column_stack([sim_color, sim_mag])
        sim_data_clean = sim_data[np.isfinite(sim_data).all(axis=1)]

        if len(sim_data_clean) < 10:
            return np.inf

        try:
            kde_sim = FFTKDE(kernel='gaussian', bw=self.bw, norm=1).fit(sim_data_clean)
            pdf_sim = kde_sim.evaluate(eval_points)

            pdf_sim = np.maximum(pdf_sim, 1e-10)
            pdf_sim = pdf_sim / np.sum(pdf_sim)

            mask = (pdf_obs > 1e-8) & (pdf_sim > 1e-8)

            if np.sum(mask) < 100:
                return np.inf

            js_distance = jensenshannon(pdf_obs[mask], pdf_sim[mask])

            if np.isnan(js_distance) or np.isinf(js_distance):
                return np.inf

            return js_distance

        except Exception as e:
            print(f"   KDE error: {e}")
            return np.inf
    
    def _extract_sim_data(self, data, filename):
        """Extract color and magnitude from simulated data"""
        # Updated column names based on your actual data structure
        color_candidates = ['BP_RP_obs', 'BP_RP', 'bp_rp', 'BP-RP', 'color']
        mag_candidates = ['G_obs', 'Gmag', 'g_mag', 'G_MAG', 'G', 'mag']
        
        color_col = None
        mag_col = None
        
        # Find color column
        for col in color_candidates:
            if col in data.columns:
                color_col = col
                break
        
        # Find magnitude column
        for col in mag_candidates:
            if col in data.columns:
                mag_col = col
                break
        
        if color_col is None or mag_col is None:
            print(f"   Warning: Could not find color/mag columns in {filename}")
            print(f"   Available columns: {list(data.columns)}")
            return None
        
        return data[color_col].values, data[mag_col].values
    
    def run_analysis(self):
        """Run the complete grid analysis"""
        print(f"\n🔬 Running grid analysis for run {self.run_number}, age {self.target_age:.1f} Gyr...")
        
        results = []
        total_models = len(self.grid_data)
        
        for i, (key, model_info) in enumerate(self.grid_data.items(), 1):
            age, dm, ebv = key
            data = model_info['data']
            filename = model_info['filename']
            filepath = model_info['filepath']
            
            # Progress indicator
            if i % 5 == 0 or i == total_models:
                print(f"   Processing model {i}/{total_models}: {filename}")
            
            try:
                # Extract simulated data
                sim_data = self._extract_sim_data(data, filename)
                if sim_data is None:
                    continue
                
                sim_color, sim_mag = sim_data
                
                # Filter to observation range
                mask = ((sim_color >= self.color_min) & (sim_color <= self.color_max) & 
                       (sim_mag >= self.mag_min) & (sim_mag <= self.mag_max))
                
                sim_color_filt = sim_color[mask]
                sim_mag_filt = sim_mag[mask]
                
                if len(sim_color_filt) < 100:
                    print(f"   Warning: Too few stars in {filename}")
                    continue
                
                # Compute distance
                js_distance = self._compute_kde_distance(sim_color_filt, sim_mag_filt)
                
                results.append({
                    'run_number': self.run_number,
                    'age': age,
                    'dm': dm,
                    'ebv': ebv,
                    'js_distance': js_distance,
                    'n_stars': len(sim_color_filt),
                    'filename': filename,
                    'filepath': filepath
                })
                
            except Exception as e:
                print(f"   Error processing {filename}: {e}")
                continue
        
        # Convert to DataFrame and sort - handle empty results
        if not results:
            print("\n❌ No valid results! Check your data columns and filtering.")
            return pd.DataFrame()
        
        results_df = pd.DataFrame(results)
        results_df = results_df.sort_values('js_distance').reset_index(drop=True)
        
        print(f"\n✅ Analysis complete! Processed {len(results_df)} models successfully")
        
        return results_df
    
    def print_best_matches(self, results_df, top_n=10):
        """Print the best matching models"""
        if len(results_df) == 0:
            print("\n❌ No results to display!")
            return pd.DataFrame()
            
        print(f"\n🏆 Top {top_n} Best Matching Models for Run {self.run_number}, Age {self.target_age:.1f} Gyr:")
        print("=" * 100)
        
        best = results_df.head(top_n)
        
        for rank, (_, row) in enumerate(best.iterrows(), 1):
            print(f"Rank {rank:2d}: Age={row['age']:4.1f} Gyr, "
                  f"DM={row['dm']:5.2f}, E(B-V)={row['ebv']:5.3f}, "
                  f"JS={row['js_distance']:.6f}, N={row['n_stars']:4d}")
        
        print("=" * 100)
        
        return best
    
    def save_results(self, results_df):
        """Save results to CSV in OUTPUT_KDE directory with age included"""
        if len(results_df) == 0:
            print("\n❌ No results to save!")
            return None
            
        output_file = f"kde_analysis_results_run_{self.run_number}_{self.target_age:.1f}.csv"
        output_path = os.path.join(self.output_path, output_file)
        
        results_df.to_csv(output_path, index=False)
        print(f"\n💾 Results saved to: {output_path}")
        
        return output_path
    
    def delete_non_best_files(self, results_df, top_n_to_keep=3):
        """Delete all CSV files except the best matching one for the target age"""
        if len(results_df) == 0:
            print("\n❌ No results to process for file deletion!")
            return

        # Get the top N matches
        top_matches = results_df.head(top_n_to_keep)
        best_filepaths = top_matches['filepath'].tolist()

        print(f"\n🗑️  Cleaning up files for run {self.run_number}, age {self.target_age:.1f} Gyr...")
        print(f"   Keeping the top {len(best_filepaths)} file(s):")
        for i, row in top_matches.iterrows():
            print(f"   - Rank {i+1}: {os.path.basename(row['filepath'])} (JS: {row['js_distance']:.6f})")

        # Get all CSV files for this specific age
        age_pattern = f"NGC6791_age{self.target_age:.1f}Gyr_*.csv"
        all_csv_files = glob.glob(f"{self.grid_path}/{age_pattern}")
        
        deleted_count = 0
        kept_count = 0
        
        for file_path in all_csv_files:
            if file_path in best_filepaths:
                kept_count += 1
            else:
                try:
                    os.remove(file_path)
                    deleted_count += 1
                except Exception as e:
                    print(f"   ❌ Error deleting {os.path.basename(file_path)}: {e}")
        
        print(f"\n📊 File cleanup summary:")
        print(f"   Files kept: {kept_count}")
        print(f"   Files deleted: {deleted_count}")

def main():
    """Main function to run the analysis"""
    # Check command line arguments
    if len(sys.argv) != 3:
        print("Usage: python kde_code_for_grid_best_match.py <run_number> <age_gyr>")
        print("Example: python kde_code_for_grid_best_match.py 20000 8.0")
        sys.exit(1)
    
    # Get run number and age from command line
    run_number = sys.argv[1]
    target_age = float(sys.argv[2])
    
    print(f"🚀 Starting NGC 6791 KDE Grid Analysis for Run {run_number}, Age {target_age:.1f} Gyr")
    print("=" * 70)
    
    try:
        # Initialize analyzer with age parameter
        analyzer = GridCMDAnalyzer(run_number, target_age=target_age)
        
        # Run analysis
        results = analyzer.run_analysis()
        
        # Only proceed if we have results
        if len(results) > 0:
            # Show results
            best_matches = analyzer.print_best_matches(results, top_n=10)
            
            # Save results
            analyzer.save_results(results)
            
            # Delete non-best files
            analyzer.delete_non_best_files(results, top_n_to_keep=4)
            
            print(f"\n🎉 Analysis complete for run {run_number}, age {target_age:.1f} Gyr!")
            
        else:
            print("\n❌ No valid results to display. Check your data format.")
            best_matches = pd.DataFrame()
        
        return analyzer, results, best_matches
        
    except Exception as e:
        print(f"\n❌ Error during analysis: {e}")
        import traceback
        traceback.print_exc()
        return None, pd.DataFrame(), pd.DataFrame()

if __name__ == "__main__":
    analyzer, results, best_matches = main()