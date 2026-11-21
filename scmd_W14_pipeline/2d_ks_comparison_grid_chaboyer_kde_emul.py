import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import jensenshannon
from scipy.interpolate import LinearNDInterpolator
from scipy.optimize import minimize
from KDEpy import FFTKDE  # Add this import
import os
import glob
import re
import sys
import argparse
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend

class CMD_2DKS_Comparison:
    def __init__(self, observed_file, simulated_file):
        """
        Initialize the 2D KS test comparison between observed and simulated CMD data.
        
        Parameters:
        -----------
        observed_file : str
            Path to the observed Gaia catalog CSV file
        simulated_file : str
            Path to the simulated CMD CSV file
        """
        self.observed_file = observed_file
        self.simulated_file = simulated_file
        self.observed_data = None
        self.simulated_data = None
        # KDE parameters - use numeric bandwidth like your working code
        self.grid_size = 1000  # Increased from 50 to match working code
        self.epsilon = 1e-6    # Decreased from 0.1 to match working code
        self.bw = 0.15         # Use fixed numeric bandwidth like working code
        
    def load_observed_data(self, color_min=1.0, color_max=1.5, mag_min=15, mag_max=19):
        """
        Load and filter the observed Gaia data.
        """
        # Load the observed data
        data = pd.read_csv(self.observed_file)
        print(f"✅ Loaded {len(data)} observed stars")
        
        # Add BP-RP color if not present
        if "BP_RP" not in data.columns:
            if "BPmag" in data.columns and "RPmag" in data.columns:
                data["BP_RP"] = data["BPmag"] - data["RPmag"]
            else:
                raise ValueError("Cannot compute BP-RP: BPmag and/or RPmag columns missing")
        
        # Filter the data based on color and magnitude
        self.observed_data = data[
            (data["BP_RP"] >= color_min) & 
            (data["BP_RP"] <= color_max) & 
            (data["Gmag"] >= mag_min) & 
            (data["Gmag"] <= mag_max)
        ].copy()
        
        print(f"✅ Filtered observed data: {len(self.observed_data)} stars")
        return self.observed_data
    
    def load_simulated_data(self, color_min=1.0, color_max=1.5, mag_min=15, mag_max=19):
        """
        Load and filter the simulated CMD data.
        """
        # Load the simulated data
        data = pd.read_csv(self.simulated_file)
        print(f"✅ Loaded {len(data)} simulated stars")
        
        # Check if columns exist, if not create them from the raw data
        if "G_obs" not in data.columns:
            data["G_obs"] = data.iloc[:, 0]  # First column is G_obs
        if "BP_obs" not in data.columns:
            data["BP_obs"] = data.iloc[:, 1]  # Second column is BP_obs
        if "RP_obs" not in data.columns:
            data["RP_obs"] = data.iloc[:, 2]  # Third column is RP_obs
        if "BP_RP_obs" not in data.columns:
            data["BP_RP_obs"] = data.iloc[:, 3]  # Fourth column is BP_RP_obs
        
        # Filter the data based on color and magnitude
        self.simulated_data = data[
            (data["BP_RP_obs"] >= color_min) & 
            (data["BP_RP_obs"] <= color_max) & 
            (data["G_obs"] >= mag_min) & 
            (data["G_obs"] <= mag_max)
        ].copy()
        
        print(f"✅ Filtered simulated data: {len(self.simulated_data)} stars")
        return self.simulated_data
    
    def _rank2d(self, points, mask=None):
        """
        Calculate 2D ranks using divide-and-conquer method.
        Adapted from your existing rankn method.
        """
        N = points.shape[0]
        N2 = N//2
        if N == 1:
            return np.array([0])
        else:
            idx = np.argpartition(points[:,0], N2)
            idxA_ = idx[:N2]
            idxA = np.zeros(N, dtype=bool)
            idxA[idxA_] = True
            
            if mask is not None:
                NAm = np.sum(idxA & mask)
                points_reduced = np.vstack((points[idxA & mask], points[~idxA & ~mask]))
            else:
                NAm = np.sum(idxA)
                points_reduced = np.vstack((points[idxA], points[~idxA]))
                
            count_points = np.zeros(points_reduced.shape[0], dtype=bool)
            count_points[:NAm] = True
            idxY = np.argsort(points_reduced[:,1])
            idxYr = np.zeros_like(idxY)
            idxYr[idxY] = np.arange(idxY.shape[0])
            count_points = count_points[idxY]
            numA = np.cumsum(count_points)[idxYr]
            
            rank = np.zeros(N, dtype=int)
            if mask is not None:
                rank[idxA] = self._rank2d(points[idxA], mask[idxA])
                rank[~idxA] = self._rank2d(points[~idxA], mask[~idxA])
                rank[~idxA & ~mask] += numA[NAm:]
            else:
                rank[idxA] = self._rank2d(points[idxA])
                rank[~idxA] = self._rank2d(points[~idxA])
                rank[~idxA] += numA[NAm:]
            return rank
    
    def calculate_2d_ks_statistic(self):
        """
        Calculate the 2D Kolmogorov-Smirnov test statistic between observed and simulated data.
        """
        if self.observed_data is None or self.simulated_data is None:
            raise ValueError("Both observed and simulated data must be loaded first")
        
        # Prepare the data arrays for comparison
        obs_points = np.column_stack([self.observed_data["BP_RP"], self.observed_data["Gmag"]])
        sim_points = np.column_stack([self.simulated_data["BP_RP_obs"], self.simulated_data["G_obs"]])
        
        # Calculate empirical CDFs using the rank method
        obs_ranks = self._rank2d(obs_points)
        obs_cdf = obs_ranks / len(obs_points)
        
        # Create interpolator for simulated data CDF
        sim_ranks = self._rank2d(sim_points)
        sim_cdf_values = sim_ranks / len(sim_points)
        sim_interpolator = LinearNDInterpolator(sim_points, sim_cdf_values, fill_value=0)
        
        # Evaluate simulated CDF at observed data points
        sim_cdf_at_obs = sim_interpolator(obs_points)
        
        # Calculate the maximum difference (KS statistic)
        differences = np.abs(obs_cdf - sim_cdf_at_obs)
        ks_statistic = np.percentile(differences, 95)

        return ks_statistic, differences
    
    def calculate_jensen_shannon_divergence(self, use_kde=True, n_bins=50):
        """
        Calculate Jensen-Shannon divergence as an alternative distance metric.
        
        Parameters:
        -----------
        use_kde : bool
            If True, use KDE method. If False, use histogram method.
        n_bins : int
            Number of bins for histogram method (ignored if use_kde=True)
        """
        if self.observed_data is None or self.simulated_data is None:
            raise ValueError("Both observed and simulated data must be loaded first")
        
        if use_kde:
            return self._compute_kde_js_divergence()
        else:
            return self._compute_histogram_js_divergence(n_bins)
    
    def _compute_kde_js_divergence(self):
        """Compute Jensen-Shannon divergence using KDE"""
        
        # Prepare data
        obs_color = self.observed_data["BP_RP"].values
        obs_mag = self.observed_data["Gmag"].values
        sim_color = self.simulated_data["BP_RP_obs"].values
        sim_mag = self.simulated_data["G_obs"].values
        
        # Determine data ranges (same as working code)
        color_min = min(obs_color.min(), sim_color.min())
        color_max = max(obs_color.max(), sim_color.max())
        mag_min = min(obs_mag.min(), sim_mag.min())
        mag_max = max(obs_mag.max(), sim_mag.max())
        
        # Define evaluation grid (use the same approach as working code)
        color_grid = np.linspace(color_min - self.epsilon, 
                                color_max + self.epsilon, 
                                self.grid_size)
        mag_grid = np.linspace(mag_min - self.epsilon, 
                              mag_max + self.epsilon, 
                              self.grid_size)
        
        # Create meshgrid and flatten
        color_mesh, mag_mesh = np.meshgrid(color_grid, mag_grid)
        eval_points = np.column_stack([color_mesh.ravel(), mag_mesh.ravel()])
        
        # Sort the evaluation points to ensure they're properly ordered
        sort_indices = np.lexsort((eval_points[:, 1], eval_points[:, 0]))
        eval_points = eval_points[sort_indices]
        
        # Prepare data for KDE (ensure it's properly formatted like working code)
        obs_data = np.column_stack([obs_color, obs_mag])
        sim_data = np.column_stack([sim_color, sim_mag])
        
        # Remove any NaN or infinite values
        obs_mask = np.isfinite(obs_data).all(axis=1)
        sim_mask = np.isfinite(sim_data).all(axis=1)
        
        obs_data_clean = obs_data[obs_mask]
        sim_data_clean = sim_data[sim_mask]
        
        if len(obs_data_clean) < 10 or len(sim_data_clean) < 10:
            return np.inf
        
        try:
            # Compute KDEs using numeric bandwidth (same as working code)
            kde_obs = FFTKDE(kernel='gaussian', bw=self.bw, norm=1).fit(obs_data_clean)
            kde_sim = FFTKDE(kernel='gaussian', bw=self.bw, norm=1).fit(sim_data_clean)
            
            # Evaluate on grid
            pdf_obs = kde_obs.evaluate(eval_points)
            pdf_sim = kde_sim.evaluate(eval_points)
            
            # Ensure PDFs are positive and normalized
            pdf_obs = np.maximum(pdf_obs, 1e-10)
            pdf_sim = np.maximum(pdf_sim, 1e-10)
            
            # Normalize to sum to 1 (create proper probability distributions)
            pdf_obs = pdf_obs / np.sum(pdf_obs)
            pdf_sim = pdf_sim / np.sum(pdf_sim)
            
            # Only use points where both PDFs have reasonable values
            mask = (pdf_obs > 1e-8) & (pdf_sim > 1e-8)
            
            if np.sum(mask) < 100:  # Need sufficient overlap
                return np.inf
            
            # Compute Jensen-Shannon distance
            js_distance = jensenshannon(pdf_obs[mask], pdf_sim[mask])
            
            # Check for valid result
            if np.isnan(js_distance) or np.isinf(js_distance):
                return np.inf
                
            return js_distance
            
        except Exception as e:
            print(f"   KDE error: {e}")
            return np.inf
    
    def _compute_histogram_js_divergence(self, n_bins=50):
        """Compute Jensen-Shannon divergence using histogram method (original)"""
        
        # Create 2D histograms
        obs_color = self.observed_data["BP_RP"]
        obs_mag = self.observed_data["Gmag"]
        sim_color = self.simulated_data["BP_RP_obs"]
        sim_mag = self.simulated_data["G_obs"]
        
        # Define common bins
        color_range = [min(obs_color.min(), sim_color.min()), 
                      max(obs_color.max(), sim_color.max())]
        mag_range = [min(obs_mag.min(), sim_mag.min()), 
                    max(obs_mag.max(), sim_mag.max())]
        
        color_bins = np.linspace(color_range[0], color_range[1], n_bins)
        mag_bins = np.linspace(mag_range[0], mag_range[1], n_bins)
        
        # Create histograms
        obs_hist, _, _ = np.histogram2d(obs_color, obs_mag, bins=[color_bins, mag_bins])
        sim_hist, _, _ = np.histogram2d(sim_color, sim_mag, bins=[color_bins, mag_bins])
        
        # Normalize to probability distributions
        obs_hist = obs_hist.flatten()
        sim_hist = sim_hist.flatten()
        obs_hist = obs_hist / obs_hist.sum()
        sim_hist = sim_hist / sim_hist.sum()
        
        # Calculate Jensen-Shannon divergence
        js_divergence = jensenshannon(obs_hist, sim_hist)
        
        return js_divergence
    
    def plot_comparison(self, color_min=1.0, color_max=1.5, mag_min=15, mag_max=19, save_path=None):
        """
        Plot observed and simulated data side by side for comparison.
        """
        if self.observed_data is None or self.simulated_data is None:
            raise ValueError("Both observed and simulated data must be loaded first")
        
        # Close any existing plots to prevent memory issues
        plt.close('all')
        
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Plot observed data (left side)
        ax1.scatter(self.observed_data["BP_RP"], self.observed_data["Gmag"], 
                   s=1, color="black", alpha=0.6, label="Observed")
        ax1.invert_yaxis()
        ax1.set_xlabel("BP - RP (mag)")
        ax1.set_ylabel("G (apparent mag)")
        ax1.set_title(f"Observed CMD (N={len(self.observed_data)})")
        ax1.set_xlim(color_min, color_max)
        ax1.set_ylim(mag_max, mag_min)
        ax1.grid(True, linestyle="--", alpha=0.3)
        
        # Plot simulated data (right side)
        ax2.scatter(self.simulated_data["BP_RP_obs"], self.simulated_data["G_obs"], 
                   s=1, color="red", alpha=0.6, label="Simulated")
        ax2.invert_yaxis()
        ax2.set_xlabel("BP - RP (mag)")
        ax2.set_ylabel("G (apparent mag)")
        ax2.set_title(f"Simulated CMD (N={len(self.simulated_data)})")
        ax2.set_xlim(color_min, color_max)
        ax2.set_ylim(mag_max, mag_min)
        ax2.grid(True, linestyle="--", alpha=0.3)
        
        # Add overall title with filename info
        if hasattr(self, 'simulated_file'):
            filename = os.path.basename(self.simulated_file)
            age, dm, ebv = parse_filename(filename)
            if age is not None:
                fig.suptitle(f"Age: {age:.1f} Gyr, DM: {dm:.2f}, E(B-V): {ebv:.2f}", fontsize=14)
        
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=120, bbox_inches='tight')
            plt.close()  # Close the figure after saving to free memory
        else:
            plt.show()
    
    def run_comparison(self, color_min=1.0, color_max=1.5, mag_min=15, mag_max=19, 
                      plot=False, save_plot=None):
        """
        Run the complete comparison analysis.
        """
        print("🔍 Starting 2D KS test comparison...")
        
        # Load data
        self.load_observed_data(color_min, color_max, mag_min, mag_max)
        self.load_simulated_data(color_min, color_max, mag_min, mag_max)
        
        # Calculate metrics
        ks_stat, differences = self.calculate_2d_ks_statistic()
        js_div = self.calculate_jensen_shannon_divergence()
        
        # Print results
        print(f"\n📊 Results:")
        print(f"   2D KS Statistic: {ks_stat:.6f}")
        print(f"   Jensen-Shannon Divergence: {js_div:.6f}")
        print(f"   Mean absolute difference: {np.mean(differences):.6f}")
        print(f"   95th percentile difference: {np.percentile(differences, 95):.6f}")
        
        # Plot if requested
        if plot:
            self.plot_comparison(color_min, color_max, mag_min, mag_max, save_plot)
        
        return {
            'ks_statistic': ks_stat,
            'js_divergence': js_div,
            'mean_diff': np.mean(differences),
            'p95_diff': np.percentile(differences, 95),
            'n_observed': len(self.observed_data),
            'n_simulated': len(self.simulated_data)
        }

def parse_filename(filename):
    """
    Parse the filename to extract age, distance modulus, and extinction.
    Example: NGC6791_age8.4Gyr_DM13.31_EBV0.13.csv
    """
    match = re.search(r'age([\d.]+)Gyr_DM([\d.]+)_EBV([\d.]+)', filename)
    if match:
        age = float(match.group(1))
        dm = float(match.group(2))
        ebv_str = match.group(3)
        # Remove any trailing periods or other non-numeric characters
        ebv_str = ebv_str.rstrip('.')
        ebv = float(ebv_str)
        return age, dm, ebv
    return None, None, None

def process_single_simulation(sim_id):
    """
    Process a single simulation directory and save results.
    """
    # File paths 
    #obs_path = "/dartfs-hpc/rc/home/8/f007058/cleaned_ngc6791_no_blue_stragglers.csv"
    obs_path = "/dartfs-hpc/rc/home/8/f007058/NGC6791_referee_isochrones_sCMDs/chaboyer_DM13.280_EBV0.200_age8.3Gyr_001.csv"
    #sim_base_path = "/dartfs-hpc/rc/lab/C/ChaboyerB/sCMD_July_11"
    sim_base_path = "/dartfs-hpc/scratch/f007058/sCMD_W14"

    #/dartfs-hpc/scratch/f007058/sCMD_W
    output_path = "/dartfs-hpc/scratch/f007058/2d_ks_W14"
    #plots_path = "/dartfs-hpc/rc/lab/C/ChaboyerB/plots_for_2d"
    plots_path = "/dartfs-hpc/scratch/f007058/plots_for_2d_W14"

    #/dartfs-hpc/scratch/f007058/plots_for_2d_W
    # Create output directories if they don't exist
    os.makedirs(output_path, exist_ok=True)
    os.makedirs(plots_path, exist_ok=True)
    
    # Process the specific simulation directory
    sim_dir = os.path.join(sim_base_path, str(sim_id))
    
    # Check if directory exists
    if not os.path.exists(sim_dir):
        print(f"❌ Directory {sim_id} does not exist at {sim_dir}")
        return None
        
    print(f"🔍 Processing directory: {sim_id}")
    print(f"   Looking in: {sim_dir}")
    
    # Get all NGC6791 files in this directory
    ngc_files = glob.glob(os.path.join(sim_dir, "NGC6791_*.csv"))
    
    if not ngc_files:
        print(f"❌ No NGC6791 files found in {sim_dir}")
        return None
    
    print(f"   Found {len(ngc_files)} NGC6791 files")
    
    results = []
    
    for ngc_file in ngc_files:
        filename = os.path.basename(ngc_file)
        print(f"  📁 Processing file: {filename}")
        
        # Parse filename to extract parameters
        age, dm, ebv = parse_filename(filename)
        if age is None:
            print(f"  ❌ Could not parse filename: {filename}")
            continue
        
        try:
            # Initialize comparison
            comparison = CMD_2DKS_Comparison(obs_path, ngc_file)
            
            # Run the analysis
            result = comparison.run_comparison(
                color_min=0.9, 
                color_max=1.7, 
                mag_min=16.25, 
                mag_max=18.0
            )
            
            # Create plot filename
            plot_filename = f"{sim_id}_{filename.replace('.csv', '')}_comparison.png"
            plot_path = os.path.join(plots_path, plot_filename)
            
            # Generate and save plot
            comparison.plot_comparison(
                color_min=0.9, 
                color_max=1.7, 
                mag_min=15.5, 
                mag_max=18.5,
                save_path=plot_path
            )
            
            # Add metadata to results
            result['simulation_id'] = str(sim_id)
            result['filename'] = filename
            result['age_gyr'] = age
            result['distance_modulus'] = dm
            result['extinction_ebv'] = ebv
            result['file_path'] = ngc_file
            result['plot_path'] = plot_path
            
            results.append(result)
            
            print(f"  ✅ KS Statistic: {result['ks_statistic']:.6f}, JS Divergence: {result['js_divergence']:.6f}")
            print(f"  📊 Plot saved: {plot_filename}")
            
        except Exception as e:
            print(f"  ❌ Error processing {filename}: {str(e)}")
            continue
    
    # Convert to DataFrame and save
    if results:
        df_results = pd.DataFrame(results)
        
        # Save results for this simulation ID
        output_file = os.path.join(output_path, f"2d_ks_results_{sim_id}.csv")
        df_results.to_csv(output_file, index=False)
        print(f"\n✅ Saved results to: {output_file}")
        
        # Find and save best fit for this simulation
        best_ks = df_results.loc[df_results['ks_statistic'].idxmin()]
        best_js = df_results.loc[df_results['js_divergence'].idxmin()]
        
        # Save best fits
        best_fits = pd.DataFrame([best_ks, best_js])
        best_fits['metric'] = ['best_ks', 'best_js']
        best_fits_file = os.path.join(output_path, f"best_fit_{sim_id}.csv")
        best_fits.to_csv(best_fits_file, index=False)
        print(f"✅ Saved best fit results to: {best_fits_file}")
        
        # Print summary
        print(f"\n📊 Summary for simulation {sim_id}:")
        print(f"   Total files processed: {len(results)}")
        print(f"   Best KS fit: {best_ks['filename']} (KS={best_ks['ks_statistic']:.6f})")
        print(f"   Best JS fit: {best_js['filename']} (JS={best_js['js_divergence']:.6f})")
        print(f"   Results saved to: {output_file}")
        print(f"   Plots saved to: {plots_path}")
        
        return df_results
    else:
        print(f"❌ No results to save for simulation {sim_id}")
        return None

def main():
    """
    Main function to handle command line arguments and run the analysis.
    """
    parser = argparse.ArgumentParser(description='Process a single simulation directory for 2D KS comparison')
    parser.add_argument('sim_id', type=int, help='Simulation ID (e.g., 20001)')
    
    args = parser.parse_args()
    
    print(f"🚀 Starting 2D KS comparison for simulation ID: {args.sim_id}")
    
    # Process the simulation
    results_df = process_single_simulation(args.sim_id)
    
    if results_df is not None:
        print(f"\n🎉 Processing complete for simulation {args.sim_id}!")
        print(f"   Check /dartfs-hpc/rc/home/8/f007058/2d_ks/ for CSV results")
        print(f"   Check /dartfs-hpc/rc/home/8/f007058/plots_for_2d/ for plots")
    else:
        print(f"\n❌ Processing failed for simulation {args.sim_id}")
        sys.exit(1)

# Main execution
if __name__ == "__main__":
    main()
