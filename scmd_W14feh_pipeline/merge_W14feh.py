import pandas as pd
import os
import glob
import numpy as np
from datetime import datetime

def combine_2d_ks_results():
    """
    Combine all individual 2D KS result CSV files into one master file.
    """
    # Paths
    results_path = "/dartfs-hpc/scratch/f007058/2d_ks_W14feh"
    output_path = "/dartfs-hpc/scratch/f007058/analysis_W14feh"
    
    print("🔍 Combining 2D KS results from individual CSV files...")
    print(f"   Looking in: {results_path}")
    
    # Find all result files
    result_files = glob.glob(os.path.join(results_path, "2d_ks_results_*.csv"))
    best_fit_files = glob.glob(os.path.join(results_path, "best_fit_*.csv"))
    
    print(f"   Found {len(result_files)} result files")
    print(f"   Found {len(best_fit_files)} best fit files")
    
    if not result_files:
        print("❌ No result files found!")
        return None
    
    # Combine all results
    all_results = []
    processed_sims = []
    failed_sims = []
    
    for result_file in result_files:
        try:
            # Extract simulation ID from filename
            filename = os.path.basename(result_file)
            sim_id = filename.replace("2d_ks_results_", "").replace(".csv", "")
            
            # Load the data
            df = pd.read_csv(result_file)
            
            if len(df) > 0:
                all_results.append(df)
                processed_sims.append(sim_id)
                print(f"  ✅ Loaded {len(df)} rows from simulation {sim_id}")
            else:
                print(f"  ⚠️  Empty file for simulation {sim_id}")
                failed_sims.append(sim_id)
                
        except Exception as e:
            print(f"  ❌ Error loading {result_file}: {str(e)}")
            failed_sims.append(sim_id)
            continue
    
    # Combine all DataFrames
    if all_results:
        combined_df = pd.concat(all_results, ignore_index=True)
        
        # Sort by simulation ID and then by KS statistic
        combined_df = combined_df.sort_values(['simulation_id', 'ks_statistic'])
        
        # Save combined results
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = os.path.join(output_path, f"combined_2d_ks_results_{timestamp}.csv")
        combined_df.to_csv(output_file, index=False)
        print(f"\n✅ Combined results saved to: {output_file}")
        
        # Create summary statistics
        create_summary_stats(combined_df, output_path, timestamp)
        
        # Combine best fits
        combine_best_fits(best_fit_files, output_path, timestamp)
        
        # Print summary
        print(f"\n📊 Summary:")
        print(f"   Total simulations processed: {len(processed_sims)}")
        print(f"   Total files analyzed: {len(combined_df)}")
        print(f"   Failed simulations: {len(failed_sims)}")
        print(f"   Unique simulation IDs: {combined_df['simulation_id'].nunique()}")
        
        if failed_sims:
            print(f"   Failed simulation IDs: {sorted(failed_sims)}")
        
        # Show some statistics
        print(f"\n📈 Quick Statistics:")
        print(f"   KS Statistic range: {combined_df['ks_statistic'].min():.6f} - {combined_df['ks_statistic'].max():.6f}")
        print(f"   JS Divergence range: {combined_df['js_divergence'].min():.6f} - {combined_df['js_divergence'].max():.6f}")
        print(f"   Age range: {combined_df['age_gyr'].min():.1f} - {combined_df['age_gyr'].max():.1f} Gyr")
        print(f"   Distance modulus range: {combined_df['distance_modulus'].min():.2f} - {combined_df['distance_modulus'].max():.2f}")
        print(f"   Extinction range: {combined_df['extinction_ebv'].min():.3f} - {combined_df['extinction_ebv'].max():.3f}")
        
        return combined_df
    else:
        print("❌ No valid results to combine!")
        return None

def create_summary_stats(df, output_path, timestamp):
    """
    Create summary statistics for the combined results.
    """
    print("📊 Creating summary statistics...")
    
    # Overall statistics
    summary_stats = {
        'total_files': len(df),
        'unique_simulations': df['simulation_id'].nunique(),
        'ks_mean': df['ks_statistic'].mean(),
        'ks_std': df['ks_statistic'].std(),
        'ks_median': df['ks_statistic'].median(),
        'ks_min': df['ks_statistic'].min(),
        'ks_max': df['ks_statistic'].max(),
        'js_mean': df['js_divergence'].mean(),
        'js_std': df['js_divergence'].std(),
        'js_median': df['js_divergence'].median(),
        'js_min': df['js_divergence'].min(),
        'js_max': df['js_divergence'].max(),
        'age_mean': df['age_gyr'].mean(),
        'age_std': df['age_gyr'].std(),
        'age_range': f"{df['age_gyr'].min():.1f} - {df['age_gyr'].max():.1f}",
        'dm_mean': df['distance_modulus'].mean(),
        'dm_std': df['distance_modulus'].std(),
        'dm_range': f"{df['distance_modulus'].min():.2f} - {df['distance_modulus'].max():.2f}",
        'ebv_mean': df['extinction_ebv'].mean(),
        'ebv_std': df['extinction_ebv'].std(),
        'ebv_range': f"{df['extinction_ebv'].min():.3f} - {df['extinction_ebv'].max():.3f}"
    }
    
    # Save summary statistics
    summary_df = pd.DataFrame([summary_stats])
    summary_file = os.path.join(output_path, f"summary_statistics_{timestamp}.csv")
    summary_df.to_csv(summary_file, index=False)
    print(f"   Summary statistics saved to: {summary_file}")
    
    # Statistics by simulation ID
    sim_stats = df.groupby('simulation_id').agg({
        'ks_statistic': ['count', 'mean', 'std', 'min', 'max'],
        'js_divergence': ['mean', 'std', 'min', 'max'],
        'age_gyr': ['mean', 'min', 'max'],
        'distance_modulus': ['mean', 'min', 'max'],
        'extinction_ebv': ['mean', 'min', 'max']
    }).round(6)
    
    # Flatten column names
    sim_stats.columns = ['_'.join(col).strip() for col in sim_stats.columns.values]
    sim_stats = sim_stats.reset_index()
    
    # Save per-simulation statistics
    sim_stats_file = os.path.join(output_path, f"per_simulation_stats_{timestamp}.csv")
    sim_stats.to_csv(sim_stats_file, index=False)
    print(f"   Per-simulation statistics saved to: {sim_stats_file}")

def combine_best_fits(best_fit_files, output_path, timestamp):
    """
    Combine all best fit files into one master file.
    """
    print("🏆 Combining best fit results...")
    
    if not best_fit_files:
        print("   No best fit files found!")
        return
    
    all_best_fits = []
    
    for best_fit_file in best_fit_files:
        try:
            # Extract simulation ID from filename
            filename = os.path.basename(best_fit_file)
            sim_id = filename.replace("best_fit_", "").replace(".csv", "")
            
            # Load the data
            df = pd.read_csv(best_fit_file)
            
            if len(df) > 0:
                all_best_fits.append(df)
                print(f"  ✅ Loaded best fits from simulation {sim_id}")
            
        except Exception as e:
            print(f"  ❌ Error loading {best_fit_file}: {str(e)}")
            continue
    
    if all_best_fits:
        combined_best_fits = pd.concat(all_best_fits, ignore_index=True)
        
        # Sort by metric and then by statistic value
        combined_best_fits = combined_best_fits.sort_values(['metric', 'ks_statistic'])
        
        # Save combined best fits
        best_fits_file = os.path.join(output_path, f"combined_best_fits_{timestamp}.csv")
        combined_best_fits.to_csv(best_fits_file, index=False)
        print(f"   Combined best fits saved to: {best_fits_file}")
        
        # Find overall best fits
        overall_best_ks = combined_best_fits[combined_best_fits['metric'] == 'best_ks'].loc[
            combined_best_fits[combined_best_fits['metric'] == 'best_ks']['ks_statistic'].idxmin()]
        overall_best_js = combined_best_fits[combined_best_fits['metric'] == 'best_js'].loc[
            combined_best_fits[combined_best_fits['metric'] == 'best_js']['js_divergence'].idxmin()]
        
        print(f"\n🏆 Overall Best Fits:")
        print(f"   Best KS: Simulation {overall_best_ks['simulation_id']}, {overall_best_ks['filename']}")
        print(f"            KS = {overall_best_ks['ks_statistic']:.6f}, JS = {overall_best_ks['js_divergence']:.6f}")
        print(f"            Age = {overall_best_ks['age_gyr']:.1f} Gyr, DM = {overall_best_ks['distance_modulus']:.2f}, E(B-V) = {overall_best_ks['extinction_ebv']:.3f}")
        
        print(f"   Best JS: Simulation {overall_best_js['simulation_id']}, {overall_best_js['filename']}")
        print(f"            KS = {overall_best_js['ks_statistic']:.6f}, JS = {overall_best_js['js_divergence']:.6f}")
        print(f"            Age = {overall_best_js['age_gyr']:.1f} Gyr, DM = {overall_best_js['distance_modulus']:.2f}, E(B-V) = {overall_best_js['extinction_ebv']:.3f}")

def find_missing_simulations():
    """
    Find which simulation IDs are missing from the results.
    """
    print("\n🔍 Checking for missing simulations...")
    
    results_path = "/dartfs-hpc/scratch/f007058/2d_ks_W14feh"
    
    # Find all result files
    result_files = glob.glob(os.path.join(results_path, "2d_ks_results_*.csv"))
    
    # Extract simulation IDs from filenames
    present_sims = set()
    for result_file in result_files:
        filename = os.path.basename(result_file)
        sim_id = filename.replace("2d_ks_results_", "").replace(".csv", "")
        try:
            present_sims.add(int(sim_id))
        except ValueError:
            continue
    
    if present_sims:
        min_sim = min(present_sims)
        max_sim = max(present_sims)
        expected_sims = set(range(min_sim, max_sim + 1))
        missing_sims = expected_sims - present_sims
        
        print(f"   Simulation range: {min_sim} - {max_sim}")
        print(f"   Present simulations: {len(present_sims)}")
        print(f"   Missing simulations: {len(missing_sims)}")
        
        if missing_sims and len(missing_sims) <= 50:  # Only print if not too many
            print(f"   Missing IDs: {sorted(missing_sims)}")
        elif missing_sims:
            print(f"   First 20 missing IDs: {sorted(list(missing_sims))[:20]}")
    
    return present_sims

def main():
    """
    Main function to combine all results.
    """
    print("🚀 Starting 2D KS results combination process...")
    
    # Find missing simulations
    present_sims = find_missing_simulations()
    
    # Combine all results
    combined_df = combine_2d_ks_results()
    
    if combined_df is not None:
        print(f"\n🎉 Results combination complete!")
        print(f"   Check /dartfs-hpc/rc/home/8/f007058/ for combined files")
    else:
        print("\n❌ Results combination failed!")

if __name__ == "__main__":
    main()