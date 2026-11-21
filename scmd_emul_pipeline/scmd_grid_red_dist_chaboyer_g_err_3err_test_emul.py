#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
simulate_scmd_5Gyr_with_errors_1D.py

This script builds a synthetic CMD (sCMD) for a 5 Gyr isochrone using Gaia data to:
  A. Load precomputed 1D empirical error distributions (magnitude-only bins).
  B. Load a 5 Gyr isochrone (with bolometric corrections via fidanka), sample stellar masses
     from a power-law IMF, assign random distances from the observed distribution,
     and add photometric errors (sampled from the 1D empirical CDFs).
  C. Adjust star counts based on completeness data.
  D. Plot the final synthetic CMD and save the simulated catalog to a CSV file.GD

Author: Georgegd
Date: 2025gd
"""

import sys
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, LinearNDInterpolator
import random

import builtins as _b
# Set DEBUG=True for verbose per-step output. Default False prints one line per stage.
DEBUG = True
def dbg(*args, **kwargs):
    if DEBUG:
        _b.print(*args, **kwargs)
def stage(msg):
    _b.print(msg, flush=True)

# Module-level caches reused across the (age, DM, EBV) grid.
_ISO_RAW_CACHE = {}      # iso_path -> parsed DataFrame
_BC_INTERP_CACHE = None  # {band: LinearNDInterpolator}
_GAIA_CACHE = {}         # gaia_file -> DataFrame
_COMPLETENESS_CACHE = {} # completeness_file -> DataFrame

_ISO_COLUMNS = [
    "EEP", "isochrone_age_yr", "initial_mass", "star_mass", "star_mdot",
    "he_core_mass", "c_core_mass", "log_L", "log_LH", "log_LHe", "log_Teff",
    "log_R", "log_g", "surface_h1", "surface_he3", "surface_he4",
    "surface_c12", "surface_o16", "log_center_T", "log_center_Rho",
    "center_gamma", "center_h1", "center_he4", "center_c12",
]

def _load_iso_raw(iso_path):
    df = _ISO_RAW_CACHE.get(iso_path)
    if df is None:
        df = pd.read_csv(iso_path, sep=r'\s+', comment="#", skiprows=10)
        df.columns = _ISO_COLUMNS
        _ISO_RAW_CACHE[iso_path] = df
    return df

def _get_bc_interpolators():
    global _BC_INTERP_CACHE
    if _BC_INTERP_CACHE is None:
        from fidanka.bolometric.load import fetch_MIST_bol_table, load_bol_table
        mist_ID = "ubv(ri)c + 2mass + kepler + hipparcos + tycho + gaia"
        _, table_folder = fetch_MIST_bol_table(mist_ID)
        bol_files = [os.path.join(table_folder, f) for f in os.listdir(table_folder) if f.startswith("feh")]
        if len(bol_files) == 0:
            raise ValueError("No bolometric correction files found!")
        bol_table_dict = load_bol_table(bol_files[3])
        bol_df = list(bol_table_dict.values())[0]
        gaia_filters = ["Gaia_G_EDR3", "Gaia_BP_EDR3", "Gaia_RP_EDR3"]
        bol_df = bol_df.dropna(subset=gaia_filters)
        _BC_INTERP_CACHE = {
            band: LinearNDInterpolator(bol_df[["Teff", "logg"]].values, bol_df[band].values)
            for band in gaia_filters
        }
    return _BC_INTERP_CACHE

# ============================
# PART A: Load Precomputed Data
# ============================

def load_precomputed_cdfs_1D(input_dir="error_cdfs_1D_June_3"):
    """
    Load precomputed 1D CDFs from text files into a dictionary.
    Each file should be named like "<band>_uncertainty_bin_<mag_bin>.txt".
    """
    error_cdfs = {}
    
    # Check if directory exists
    if not os.path.exists(input_dir):
        dbg(f"⚠️ Error directory {input_dir} not found!")
        return error_cdfs
    
    # Print all files in directory for debugging
    dbg(f"Files in {input_dir}:")
    all_files = os.listdir(input_dir)
    dbg(f"Found {len(all_files)} files")
    
    # Process each file with proper error handling
    successful_files = 0
    for file_name in all_files:
        if file_name.endswith(".txt"):
            try:
                # Parse filename parts
                parts = file_name.split("_")
                band = parts[0]  # g, bp, rp
                
                # Extract magnitude bin number from filename
                # Format is like: g_uncertainty_bin_0.txt
                mag_bin = int(parts[3].split('.')[0])
                
                key = mag_bin  # Only magnitude bin, no distance bin
                
                # Load file data
                file_path = os.path.join(input_dir, file_name)
                with open(file_path, "r") as f:
                    lines = f.readlines()[2:]  # Skip header lines
                    
                    # Parse data lines
                    errors = []
                    cdfs = []
                    for line in lines:
                        if line.strip():  # Skip empty lines
                            parts = line.strip().split()
                            if len(parts) >= 2:
                                errors.append(float(parts[0]))
                                cdfs.append(float(parts[1]))
                    
                    # Store in dictionary
                    if key not in error_cdfs:
                        error_cdfs[key] = {}
                    
                    error_cdfs[key][f"sorted_errors_{band}"] = np.array(errors)
                    error_cdfs[key][f"cdf_{band}"] = np.array(cdfs)
                    successful_files += 1
                    
            except Exception as e:
                dbg(f"⚠️ Error processing file {file_name}: {str(e)}")
    
    # Print summary
    dbg(f"✅ Successfully loaded {successful_files} error files for {len(error_cdfs)} magnitude bins")
    
    # Print some sample data for verification
    if error_cdfs:
        sample_key = list(error_cdfs.keys())[0]
        dbg(f"Sample mag bin {sample_key}:")
        for key in sorted(error_cdfs[sample_key].keys()):
            data = error_cdfs[sample_key][key]
            dbg(f"  - {key}: {len(data)} points, range: {min(data):.6f} to {max(data):.6f}")
    
    return error_cdfs

def load_completeness(file_path):
    """
    Load the completeness data from a CSV file.
    """
    completeness_data = pd.read_csv(file_path)
    dbg(f"✅ Loaded completeness data from {file_path}")
    return completeness_data

# ============================
# PART B: Isochrone Processing and Synthetic CMD Simulation
# ============================

def apply_extinction(magnitudes_dict, ebv=0.16):
    """
    Apply extinction to Gaia magnitudes using standard extinction coefficients.
    
    Parameters:
    -----------
    magnitudes_dict : dict
        Dictionary with keys "Gaia_G_EDR3", "Gaia_BP_EDR3", "Gaia_RP_EDR3"
    ebv : float
        E(B-V) extinction value in magnitudes
        
    Returns:
    --------
    reddened_magnitudes : dict
        Dictionary with reddened magnitudes
    """
    # Gaia extinction coefficients from Gaia Collaboration et al. (2018)
    # A_G/E(B-V), A_BP/E(B-V), A_RP/E(B-V)
    extinction_coeffs = {
        "Gaia_G_EDR3": 2.740,   # A_G/E(B-V)
        "Gaia_BP_EDR3": 3.374,  # A_BP/E(B-V) 
        "Gaia_RP_EDR3": 2.035   # A_RP/E(B-V)
    }
    
    reddened_magnitudes = {}
    for band in magnitudes_dict:
        A_band = extinction_coeffs[band] * ebv
        reddened_magnitudes[band] = magnitudes_dict[band] + A_band
    
    return reddened_magnitudes

def load_isochrone_with_bol_corr(iso_path, target_age_years, apparent_distance_modulus=13.0, ebv=0.16):
    """
    Load isochrone for a specific age and convert to apparent magnitudes with apparent distance modulus.
    
    Parameters:
    -----------
    iso_path : str
        Path to the isochrone file
    target_age_years : float
        Target age in years (not log10!)
    apparent_distance_modulus : float
        Apparent distance modulus (includes extinction effects) from literature
    ebv : float
        E(B-V) extinction value (used only for color calculations)
    """
    iso_df = _load_iso_raw(iso_path)

    # Filter for the desired age
    iso_df = iso_df[np.isclose(iso_df["isochrone_age_yr"], target_age_years, rtol=0.01)]
    if iso_df.empty:
        raise ValueError(f"No isochrone data found for age = {target_age_years:.0e} years")

    iso_Teff = 10 ** iso_df["log_Teff"]
    iso_logg = iso_df["log_g"]
    iso_logL = iso_df["log_L"]
    iso_ages = iso_df["isochrone_age_yr"]

    bc_interpolators = _get_bc_interpolators()
    gaia_filters = ["Gaia_G_EDR3", "Gaia_BP_EDR3", "Gaia_RP_EDR3"]
    M_bol = 4.75 - 2.5 * iso_logL
    magnitudes = {}
    for band in gaia_filters:
        BC = bc_interpolators[band](iso_Teff, iso_logg)
        magnitudes[band] = M_bol - BC
    
    # ⭐ APPLY EXTINCTION TO ABSOLUTE MAGNITUDES ⭐
    reddened_magnitudes = apply_extinction(magnitudes, ebv)
    
    # ⭐ CONVERT TO APPARENT MAGNITUDES USING TRUE DISTANCE MODULUS ⭐
    # Calculate true distance modulus by removing extinction from apparent DM
    extinction_coeffs = {
        "Gaia_G_EDR3": 2.740,   # A_G/E(B-V)
        "Gaia_BP_EDR3": 3.374,  # A_BP/E(B-V) 
        "Gaia_RP_EDR3": 2.035   # A_RP/E(B-V)
    }
    
    # Calculate true distance modulus for each band
    A_G = extinction_coeffs["Gaia_G_EDR3"] * ebv
    true_distance_modulus = apparent_distance_modulus - A_G
    
    iso_df_out = pd.DataFrame(reddened_magnitudes)
    iso_df_out["initial_mass"] = iso_df["initial_mass"]
    iso_df_out["log_Teff"] = np.log10(iso_Teff)
    iso_df_out["log_L"] = iso_logL
    iso_df_out["isochrone_age_yr"] = iso_ages
    
    # Add the TRUE distance modulus (extinction already applied above)
    iso_df_out["Gaia_G_EDR3"] += true_distance_modulus
    iso_df_out["Gaia_BP_EDR3"] += true_distance_modulus  
    iso_df_out["Gaia_RP_EDR3"] += true_distance_modulus
    
    # Color calculation (extinction affects color!)
    iso_df_out["BP-RP"] = iso_df_out["Gaia_BP_EDR3"] - iso_df_out["Gaia_RP_EDR3"]
    
    age_gyr = target_age_years / 1e9
    dbg(f"🔧 Isochrone loaded for {age_gyr:.1f} Gyr ({target_age_years:.0e} years)")
    dbg(f"   Apparent distance modulus: {apparent_distance_modulus:.2f} mag")
    dbg(f"   True distance modulus: {true_distance_modulus:.2f} mag")
    dbg(f"   Applied extinction E(B-V): {ebv:.2f} mag (A_G = {A_G:.2f} mag)")
    dbg(f"   G range: {iso_df_out['Gaia_G_EDR3'].min():.2f} to {iso_df_out['Gaia_G_EDR3'].max():.2f}")
    dbg(f"   BP-RP range: {iso_df_out['BP-RP'].min():.2f} to {iso_df_out['BP-RP'].max():.2f}")
    
    return iso_df_out

def find_msto_point(iso_df_out):
    """
    Find the Main Sequence Turn-Off point in the isochrone.
    The MSTO is identified as the bluest point (minimum BP-RP) in the isochrone.
    
    Parameters:
    -----------
    iso_df_out : pandas.DataFrame
        The isochrone DataFrame with Gaia apparent magnitudes
        
    Returns:
    --------
    msto_mag : float
        The apparent G magnitude at the MSTO point
    msto_color : float
        The BP-RP color at the MSTO point
    """
    # Find the bluest point (minimum BP-RP color)
    color = iso_df_out["Gaia_BP_EDR3"] - iso_df_out["Gaia_RP_EDR3"]
    min_color_idx = np.argmin(color)
    
    msto_mag = iso_df_out["Gaia_G_EDR3"].iloc[min_color_idx]
    msto_color = color.iloc[min_color_idx]
    
    dbg(f"✅ Found MSTO point: G(apparent) = {msto_mag:.3f}, BP-RP = {msto_color:.3f}")
    
    return msto_mag, msto_color

def simulate_stars(iso_df_out, nstars=100000, IMF_slope=-1.02):
    """
    Simulate stars - everything is already in apparent magnitudes from the isochrone
    """
    mass_iso = iso_df_out["initial_mass"].values
    m_min = np.min(mass_iso)
    m_max = np.max(mass_iso)
    
    def sample_imf(n, m_min, m_max, slope):
        r = np.random.random(n)
        exponent = slope + 1.0
        return (r * (m_max ** exponent - m_min ** exponent) + m_min ** exponent) ** (1/exponent)
    
    primary_masses = sample_imf(nstars, m_min, m_max, IMF_slope)
    sort_idx = np.argsort(mass_iso)
    m_sorted = mass_iso[sort_idx]
    
    # These are now apparent magnitudes!
    G_iso = iso_df_out["Gaia_G_EDR3"].values[sort_idx]
    BP_iso = iso_df_out["Gaia_BP_EDR3"].values[sort_idx]
    RP_iso = iso_df_out["Gaia_RP_EDR3"].values[sort_idx]
    
    G_interp = interp1d(m_sorted, G_iso, kind="linear", fill_value="extrapolate")
    BP_interp = interp1d(m_sorted, BP_iso, kind="linear", fill_value="extrapolate")
    RP_interp = interp1d(m_sorted, RP_iso, kind="linear", fill_value="extrapolate")
    
    # Sample apparent magnitudes directly
    sim_G = G_interp(primary_masses)
    sim_BP = BP_interp(primary_masses)
    sim_RP = RP_interp(primary_masses)
    
    dbg(f"🔧 Simulated {nstars} stars with apparent magnitudes")
    dbg(f"   G range: {sim_G.min():.2f} to {sim_G.max():.2f}")
    
    # Apply binary stars with 2% fraction.
    binary_fraction = 0.0
    for i in range(nstars):
        if np.random.random() <= binary_fraction:
            q = 0.5 + 0.5 * np.random.random()
            secondary_mass = q * primary_masses[i]
            sec_G = G_interp(secondary_mass)
            sec_BP = BP_interp(secondary_mass)
            sec_RP = RP_interp(secondary_mass)
            combined_flux_G = 10 ** (-0.4 * sim_G[i]) + 10 ** (-0.4 * sec_G)
            sim_G[i] = -2.5 * np.log10(combined_flux_G)
            combined_flux_BP = 10 ** (-0.4 * sim_BP[i]) + 10 ** (-0.4 * sec_BP)
            sim_BP[i] = -2.5 * np.log10(combined_flux_BP)
            combined_flux_RP = 10 ** (-0.4 * sim_RP[i]) + 10 ** (-0.4 * sec_RP)
            sim_RP[i] = -2.5 * np.log10(combined_flux_RP)
    
    sim_color = sim_BP - sim_RP
    return primary_masses, sim_G, sim_BP, sim_RP, sim_color

def add_photometric_errors_1D(sim_G, sim_BP, sim_RP, error_cdfs_1D, mag_bins):
    """
    Add realistic Gaia‐like photometric errors to synthetic magnitudes using 1D magnitude-only bins.
    VECTORIZED VERSION for better performance.
    """
    import numpy as np

    dbg("🔧 Adding photometric errors to apparent magnitudes (vectorized)...")
    
    # Convert to numpy arrays
    sim_G = np.array(sim_G, dtype=float)
    sim_BP = np.array(sim_BP, dtype=float)
    sim_RP = np.array(sim_RP, dtype=float)
    
    # Prepare outputs
    sim_G_obs = sim_G.copy()
    sim_BP_obs = sim_BP.copy()
    sim_RP_obs = sim_RP.copy()

    # Use apparent G magnitudes for binning
    G_for_binning = sim_G
    Nmagbins = len(mag_bins) - 1
    
    # Vectorized bin assignment
    bin_indices = np.clip(
        np.digitize(G_for_binning, mag_bins, right=False) - 1,
        0, Nmagbins - 1
    )
    
    # Create masks for each magnitude bin that has error data
    available_bins = set(error_cdfs_1D.keys())
    success_mask = np.zeros(len(sim_G), dtype=bool)
    
    # Process each bin that has error data
    for m_bin in available_bins:
        if m_bin not in error_cdfs_1D:
            continue
            
        # Find stars in this magnitude bin
        bin_mask = (bin_indices == m_bin)
        n_stars_in_bin = np.sum(bin_mask)
        
        if n_stars_in_bin == 0:
            continue
        
        # Apply G band errors
        if "sorted_errors_g" in error_cdfs_1D[m_bin] and "cdf_g" in error_cdfs_1D[m_bin]:
            se_g = error_cdfs_1D[m_bin]["sorted_errors_g"]
            cdf_g = error_cdfs_1D[m_bin]["cdf_g"]
            if len(se_g) > 0 and len(cdf_g) > 0:
                # Sample random values for all stars in this bin
                random_vals = np.random.rand(n_stars_in_bin)
                # Interpolate to get error values
                sigma_g_vals = np.interp(random_vals, cdf_g, se_g)
                # Add Gaussian noise
                noise_g = np.random.normal(0, sigma_g_vals)
                sim_G_obs[bin_mask] += noise_g
        
        # Apply BP band errors
        if "sorted_errors_bp" in error_cdfs_1D[m_bin] and "cdf_bp" in error_cdfs_1D[m_bin]:
            se_bp = error_cdfs_1D[m_bin]["sorted_errors_bp"]
            cdf_bp = error_cdfs_1D[m_bin]["cdf_bp"]
            if len(se_bp) > 0 and len(cdf_bp) > 0:
                random_vals = np.random.rand(n_stars_in_bin)
                sigma_bp_vals = np.interp(random_vals, cdf_bp, se_bp)
                noise_bp = np.random.normal(0, sigma_bp_vals)
                sim_BP_obs[bin_mask] += noise_bp
        
        # Apply RP band errors
        if "sorted_errors_rp" in error_cdfs_1D[m_bin] and "cdf_rp" in error_cdfs_1D[m_bin]:
            se_rp = error_cdfs_1D[m_bin]["sorted_errors_rp"]
            cdf_rp = error_cdfs_1D[m_bin]["cdf_rp"]
            if len(se_rp) > 0 and len(cdf_rp) > 0:
                random_vals = np.random.rand(n_stars_in_bin)
                sigma_rp_vals = np.interp(random_vals, cdf_rp, se_rp)
                noise_rp = np.random.normal(0, sigma_rp_vals)
                sim_RP_obs[bin_mask] += noise_rp
        
        # Mark these stars as successfully processed
        success_mask[bin_mask] = True
    
    # Apply fallback errors to stars without empirical error data
    fallback_mask = ~success_mask
    n_fallback = np.sum(fallback_mask)
    
    if n_fallback > 0:
        # Vectorized fallback error calculation
        G_fallback = G_for_binning[fallback_mask]
        sigma_fb = 0.01 * (1.0 + 0.1 * (G_fallback - 13))
        
        # Add fallback noise
        sim_G_obs[fallback_mask] += np.random.normal(0, sigma_fb)
        sim_BP_obs[fallback_mask] += np.random.normal(0, sigma_fb * 1.2)
        sim_RP_obs[fallback_mask] += np.random.normal(0, sigma_fb * 1.1)
    
    sim_color_obs = sim_BP_obs - sim_RP_obs
    
    success_count = np.sum(success_mask)
    fallback_count = n_fallback
    
    dbg(f"✅ Empirical errors applied to {success_count} stars; fallback used for {fallback_count} stars.")
    dbg(f"   Final G range: {sim_G_obs.min():.2f} to {sim_G_obs.max():.2f}")
    
    return sim_G_obs, sim_BP_obs, sim_RP_obs, sim_color_obs

def add_perpendicular_errors(sim_G, sim_BP, sim_RP, iso_df_out, perp_error_stats):
    """
    Add perpendicular errors to synthetic magnitudes based on RP bins and isochrone direction.
    VECTORIZED VERSION for better performance.
    """
    dbg("🔧 Adding perpendicular errors based on RP bins (vectorized)...")
    
    # Convert to numpy arrays
    sim_G = np.array(sim_G, dtype=float)
    sim_BP = np.array(sim_BP, dtype=float) 
    sim_RP = np.array(sim_RP, dtype=float)
    
    # Calculate isochrone slopes
    slope_interpolator = calculate_isochrone_slope(iso_df_out)
    
    # Initialize output arrays
    sim_G_obs = sim_G.copy()
    sim_BP_obs = sim_BP.copy()
    sim_RP_obs = sim_RP.copy()
    
    # Get RP bin information
    rp_left = perp_error_stats['RP_left'].values
    rp_right = perp_error_stats['RP_right'].values
    std_scatter = perp_error_stats['std_scatter'].values#*1.1
    
    # Vectorized bin assignment for all stars
    n_stars = len(sim_RP)
    bin_assignments = np.full(n_stars, -1, dtype=int)  # -1 means no bin found
    
    # Find bin for each star (vectorized)
    for bin_idx in range(len(rp_left)):
        mask = (sim_RP >= rp_left[bin_idx]) & (sim_RP < rp_right[bin_idx])
        bin_assignments[mask] = bin_idx
    
    # Create mask for stars that have a valid bin assignment
    valid_mask = bin_assignments >= 0
    n_valid = np.sum(valid_mask)
    
    if n_valid > 0:
        # Get data for valid stars only
        valid_rp = sim_RP[valid_mask]
        valid_bins = bin_assignments[valid_mask]
        
        # Get scatter values for each valid star
        sigma_perp_vals = std_scatter[valid_bins]
        
        # Get isochrone slopes for all valid RP values at once
        iso_slopes = slope_interpolator(valid_rp)
        
        # Calculate perpendicular slopes (vectorized)
        # Avoid division by zero
        perp_slopes = np.where(np.abs(iso_slopes) > 1e-6, 
                              -1.0 / iso_slopes, 
                              np.sign(iso_slopes) * 1e6)
        
        # Sample perpendicular error magnitudes
        error_magnitudes = np.random.normal(0, sigma_perp_vals)
        
        # Calculate unit perpendicular vectors (vectorized)
        perp_norms = np.sqrt(1 + perp_slopes**2)
        unit_perp_rp = 1.0 / perp_norms
        unit_perp_bp_rp = perp_slopes / perp_norms;
        
        # Apply perpendicular errors (vectorized)
        rp_errors = error_magnitudes * unit_perp_rp
        bp_rp_errors = error_magnitudes * unit_perp_bp_rp;
        
        # Update magnitudes for valid stars
        sim_RP_obs[valid_mask] += rp_errors
        sim_BP_obs[valid_mask] += rp_errors + bp_rp_errors  # BP = RP + (BP-RP)
        
        # Add correlated G band errors (scaled version of perpendicular error)
        g_errors = np.random.normal(0, sigma_perp_vals * 0.8)
        sim_G_obs[valid_mask] += g_errors
    
    # Apply fallback errors to stars outside RP bin range
    fallback_mask = ~valid_mask
    n_fallback = np.sum(fallback_mask)
    
    if n_fallback > 0:
        sigma_fb = 0.05  # Default fallback error
        sim_G_obs[fallback_mask] += np.random.normal(0, sigma_fb, n_fallback)
        sim_BP_obs[fallback_mask] += np.random.normal(0, sigma_fb, n_fallback)
        sim_RP_obs[fallback_mask] += np.random.normal(0, sigma_fb, n_fallback)
    
    sim_color_obs = sim_BP_obs - sim_RP_obs
    
    dbg(f"✅ Applied perpendicular errors to {n_valid} stars")
    dbg(f"   Fallback errors for {n_fallback} stars outside RP bin range")
    dbg(f"   Final G range: {sim_G_obs.min():.2f} to {sim_G_obs.max():.2f}")
    dbg(f"   Final BP-RP range: {sim_color_obs.min():.2f} to {sim_color_obs.max():.2f}")
    
    return sim_G_obs, sim_BP_obs, sim_RP_obs, sim_color_obs

def load_merged_scatter_data(file_path):
    """
    Load the merged scatter statistics from the CSV file.
    
    Parameters:
    -----------
    file_path : str
        Path to merged_bins_analysis.csv
        
    Returns:
    --------
    scatter_stats : pandas.DataFrame
        Scatter statistics with G magnitude bins and std_scatter values
    """
    scatter_stats = pd.read_csv(file_path)
    
    dbg(f"✅ Loaded merged scatter statistics:")
    dbg(f"   Total G mag bins: {len(scatter_stats)}")
    dbg(f"   G mag range: {scatter_stats['g_mag_mean'].min():.2f} to {scatter_stats['g_mag_mean'].max():.2f}")
    dbg(f"   Scatter range: {scatter_stats['std_scatter'].min():.4f} to {scatter_stats['std_scatter'].max():.4f}")
    
    return scatter_stats

def add_g_magnitude_scatter(sim_G, sim_BP, sim_RP, iso_df_out, scatter_stats):
    """
    Add perpendicular scatter to synthetic magnitudes based on G magnitude bins and isochrone direction.
    Uses std_scatter from merged_bins_analysis.csv and injects perpendicular to the isochrone
    in G vs RP-BP color-magnitude diagram.
    Maximum perpendicular distance is limited to 0.08 in normalized coordinates.
    """
    dbg("🔧 Adding G magnitude-based scatter perpendicular to isochrone (vectorized)...")
    
    # Convert to numpy arrays
    sim_G = np.array(sim_G, dtype=float)
    sim_BP = np.array(sim_BP, dtype=float) 
    sim_RP = np.array(sim_RP, dtype=float)
    
    # Calculate isochrone slopes in G vs RP-BP space
    slope_interpolator = calculate_isochrone_slope_g_vs_color(iso_df_out)
    
    # Initialize output arrays
    sim_G_obs = sim_G.copy()
    sim_BP_obs = sim_BP.copy()
    sim_RP_obs = sim_RP.copy()
    
    # Get G magnitude bin information from scatter_stats
    g_mag_means = scatter_stats['g_mag_mean'].values
    std_scatter_vals = scatter_stats['std_scatter'].values
    bin_widths = scatter_stats['bin_width'].values
    
    # Calculate normalization scales from isochrone data
    color_iso = iso_df_out["Gaia_BP_EDR3"] - iso_df_out["Gaia_RP_EDR3"]
    g_iso = iso_df_out["Gaia_G_EDR3"]
    color_scale = color_iso.max() - color_iso.min()  # Color axis range
    mag_scale = g_iso.max() - g_iso.min()           # Magnitude axis range
    
    # Maximum allowed perpendicular distance in normalized units
    max_perp_distance_normalized = 0.05
    
    dbg(f"   Normalization scales: color_scale={color_scale:.3f}, mag_scale={mag_scale:.3f}")
    dbg(f"   Max perpendicular distance: {max_perp_distance_normalized} (normalized units)")
    
    # Vectorized bin assignment for all stars based on G magnitude
    n_stars = len(sim_G)
    bin_assignments = np.full(n_stars, -1, dtype=int)  # -1 means no bin found
    
    # Find bin for each star based on G magnitude (vectorized)
    for bin_idx in range(len(scatter_stats)):
        g_center = g_mag_means[bin_idx]
        bin_width = bin_widths[bin_idx]
        g_left = g_center - bin_width/2
        g_right = g_center + bin_width/2
        
        mask = (sim_G >= g_left) & (sim_G < g_right)
        bin_assignments[mask] = bin_idx
    
    # Only process stars with valid bins
    valid_mask = bin_assignments >= 0
    n_valid = np.sum(valid_mask)
    
    if n_valid > 0:
        # Get valid star data
        valid_g = sim_G[valid_mask]
        valid_rp_bp = sim_BP[valid_mask] - sim_RP[valid_mask]  # RP-BP color
        valid_bins = bin_assignments[valid_mask]
        
        # Get scatter values for each valid star
        sigma_scatter_vals = std_scatter_vals[valid_bins]
        
        # Get isochrone slopes for all valid G values at once
        iso_slopes = slope_interpolator(valid_g)
        
        # Calculate perpendicular slopes (vectorized)
        perp_slopes = np.where(np.abs(iso_slopes) > 1e-6, 
                              -1.0 / iso_slopes, 
                              np.sign(iso_slopes) * 1e6)
        
        # Sample perpendicular error magnitudes using std_scatter
        error_magnitudes = np.random.normal(0, sigma_scatter_vals)
        
        # ⭐ APPLY MAXIMUM DISTANCE CONSTRAINT IN NORMALIZED COORDINATES ⭐
        # Calculate unit perpendicular vectors in G vs BP-RP space (unnormalized)
        perp_norms = np.sqrt(1 + perp_slopes**2)
        unit_perp_g = 1.0 / perp_norms
        unit_perp_color = perp_slopes / perp_norms;
        
        # Calculate actual movements in physical space
        g_errors = error_magnitudes * unit_perp_g
        color_errors = error_magnitudes * unit_perp_color;
        
        # Normalize the perpendicular distances
        g_errors_norm = g_errors / mag_scale
        color_errors_norm = color_errors / color_scale
        distances_norm = np.sqrt(g_errors_norm**2 + color_errors_norm**2)
        
        # Clip to maximum normalized distance
        scale_factors = np.where(distances_norm > max_perp_distance_normalized,
                                max_perp_distance_normalized / distances_norm,
                                1.0)
        
        g_errors *= scale_factors
        color_errors *= scale_factors
        
        # Update G magnitudes
        sim_G_obs[valid_mask] += g_errors
        
        # Apply color errors to both BP and RP to change RP-BP
        bp_errors = color_errors * 0.5
        rp_errors = -color_errors * 0.5;
        
        sim_BP_obs[valid_mask] += bp_errors
        sim_RP_obs[valid_mask] += rp_errors
        
        # Count how many stars were clipped
        n_clipped = np.sum(scale_factors < 1.0)
        if n_clipped > 0:
            dbg(f"   Applied normalized distance limit: {n_clipped} stars clipped to max distance {max_perp_distance_normalized}")
    
    # Apply fallback errors to stars outside G bin range
    fallback_mask = ~valid_mask
    n_fallback = np.sum(fallback_mask)
    
    if n_fallback > 0:
        sigma_fb = 0.02  # Default fallback scatter
        fallback_errors_g = np.random.normal(0, sigma_fb, n_fallback)
        fallback_errors_bp = np.random.normal(0, sigma_fb * 0.5, n_fallback)
        fallback_errors_rp = np.random.normal(0, sigma_fb * 0.5, n_fallback)
        
        # Normalize and clip fallback errors
        fallback_errors_g = np.clip(fallback_errors_g, -max_perp_distance_normalized, max_perp_distance_normalized)
        fallback_errors_bp = np.clip(fallback_errors_bp, -max_perp_distance_normalized, max_perp_distance_normalized)
        fallback_errors_rp = np.clip(fallback_errors_rp, -max_perp_distance_normalized, max_perp_distance_normalized)
        
        sim_G_obs[fallback_mask] += fallback_errors_g
        sim_BP_obs[fallback_mask] += fallback_errors_bp
        sim_RP_obs[fallback_mask] += fallback_errors_rp
    
    sim_color_obs = sim_BP_obs - sim_RP_obs
    
    dbg(f"✅ Applied G magnitude-based scatter to {n_valid} stars")
    dbg(f"   Fallback errors for {n_fallback} stars outside G bin range")
    dbg(f"   Final G range: {sim_G_obs.min():.2f} to {sim_G_obs.max():.2f}")
    dbg(f"   Final RP-BP range: {sim_color_obs.min():.2f} to {sim_color_obs.max():.2f}")
    
    return sim_G_obs, sim_BP_obs, sim_RP_obs, sim_color_obs

def calculate_isochrone_slope_g_vs_color(iso_df_out):
    """
    Calculate the local slope of the isochrone in G vs BP-RP color space.
    
    Parameters:
    -----------
    iso_df_out : pandas.DataFrame
        Isochrone data with Gaia magnitudes
        
    Returns:
    --------
    slope_interpolator : function
        Function that returns slope as function of G magnitude
    """
    # Get G magnitudes and colors from isochrone
    G_iso = iso_df_out["Gaia_G_EDR3"].values
    color_iso = iso_df_out["Gaia_BP_EDR3"].values - iso_df_out["Gaia_RP_EDR3"].values
    
    # Sort by G magnitude
    sort_idx = np.argsort(G_iso)
    G_sorted = G_iso[sort_idx]
    color_sorted = color_iso[sort_idx]
    
    # Calculate slopes using finite differences
    slopes = []
    G_centers = []
    
    for i in range(1, len(G_sorted) - 1):
        # Use neighboring parts to calculate slope
        dG = G_sorted[i+1] - G_sorted[i-1]
        dcolor = color_sorted[i+1] - color_sorted[i-1]
        
        if abs(dG) > 1e-6:  # Avoid division by zero
            slope = dcolor / dG  # dcolor/dG
            slopes.append(slope)
            G_centers.append(G_sorted[i])
    
    # Handle edge cases
    if len(slopes) > 0:
        # Add endpoints with nearest neighbor slopes
        G_centers = [G_sorted[0]] + G_centers + [G_sorted[-1]]
        slopes = [slopes[0]] + slopes + [slopes[-1]]
    else:
        # Fallback to constant slope
        G_centers = [G_sorted[0], G_sorted[-1]]
        slopes = [0.0, 0.0]
    
    # Create interpolator
    slope_interpolator = interp1d(G_centers, slopes, kind='linear', 
                                 bounds_error=False, fill_value='extrapolate')
    
    dbg(f"✅ Created isochrone slope interpolator for G vs BP-RP")
    dbg(f"   G range: {min(G_centers):.2f} to {max(G_centers):.2f}")
    dbg(f"   Slope range: {min(slopes):.3f} to {max(slopes):.3f}")
    
    return slope_interpolator

def load_regional_scatter_data(red_scatter_file, blue_scatter_file):
    """
    Load regional scatter statistics from CSV files.
    
    Parameters:
    -----------
    red_scatter_file : str
        Path to red region bin statistics CSV
    blue_scatter_file : str
        Path to blue region bin statistics CSV
        
    Returns:
    --------
    red_stats : pandas.DataFrame
        Red region scatter statistics with columns:
        - bin_center, bin_left, bin_right, bin_width
        - n_stars, g_mean, g_median, bp_rp_mean, bp_rp_median
        - perp_mean, perp_median, perp_std, perp_mad
        - perp_16th, perp_84th, perp_scatter_68
    blue_stats : pandas.DataFrame
        Blue region scatter statistics (same structure)
    """
    # Load red region stats
    red_stats = pd.read_csv(red_scatter_file)
    blue_stats = pd.read_csv(blue_scatter_file)
    
    dbg(f"✅ Loaded regional scatter statistics:")
    dbg(f"   Red region: {len(red_stats)} G magnitude bins")
    dbg(f"     G range: {red_stats['bin_left'].min():.2f} to {red_stats['bin_right'].max():.2f}")
    dbg(f"     Stars per bin: {red_stats['n_stars'].min()}-{red_stats['n_stars'].max()}")
    dbg(f"     Perpendicular std: {red_stats['perp_std'].min():.4f} to {red_stats['perp_std'].max():.4f} mag")
    
    dbg(f"   Blue region: {len(blue_stats)} G magnitude bins")
    dbg(f"     G range: {blue_stats['bin_left'].min():.2f} to {blue_stats['bin_right'].max():.2f}")
    dbg(f"     Stars per bin: {blue_stats['n_stars'].min()}-{blue_stats['n_stars'].max()}")
    dbg(f"     Perpendicular std: {blue_stats['perp_std'].min():.4f} to {blue_stats['perp_std'].max():.4f} mag")
    
    return red_stats, blue_stats

def process_single_isochrone_grid_quota(age_years, apparent_distance_modulus, ebv, 
                                       isochrone_file, gaia_file, red_stats, blue_stats,
                                       completeness_file, output_dir="grid_results"):
    """
    Process a single combination using apparent distance modulus from literature.
    Updated to use regional G magnitude-based scatter.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    age_gyr = age_years / 1e9
    
    dbg(f"\n{'='*80}")
    dbg(f"🔧 PROCESSING COMBINATION")
    dbg(f"   Age: {age_years:.0e} years ({age_gyr:.1f} Gyr)")
    dbg(f"   Distance modulus: {apparent_distance_modulus:.2f} mag")
    dbg(f"   E(B-V): {ebv:.3f} mag")
    dbg(f"{'='*80}")
    
    try:
        # Load the isochrone with apparent distance modulus
        iso_df_out = load_isochrone_with_bol_corr(isochrone_file, age_years, apparent_distance_modulus, ebv)
        
        # Find MSTO in the apparent magnitude isochrone
        msto_mag, msto_color = find_msto_point(iso_df_out)
        
        # Generate synthetic stars using quota sampling
        dbg(f"🔧 Generating stars with quota sampling...")
        synthetic_catalog = quota_sample_stars(iso_df_out, gaia_file, completeness_file, imf_slope=-1.02)
        
        if len(synthetic_catalog) == 0:
            dbg("❌ No stars generated!")
            return False, 0, np.nan, np.nan
        
        # Add regional G magnitude-based scatter
        dbg(f"🔧 Adding regional G magnitude-based scatter...")
        sim_G_obs, sim_BP_obs, sim_RP_obs, sim_color_obs = add_regional_g_magnitude_scatter(
            synthetic_catalog['G'], synthetic_catalog['BP'], synthetic_catalog['RP'], 
            iso_df_out, red_stats, blue_stats
        )
        
        # Create final catalog with minimal columns - only observed magnitudes
        final_catalog = pd.DataFrame({
            "G_obs": sim_G_obs,
            "BP_obs": sim_BP_obs,
            "RP_obs": sim_RP_obs,
            "BP_RP_obs": sim_color_obs
        })
        
        dbg(f"📊 Final catalog: {len(final_catalog):,} stars")
        
        # Generate filename
        csv_filename = generate_filename(age_years, apparent_distance_modulus, ebv, "csv")
        
        # Save CSV
        csv_path = os.path.join(output_dir, csv_filename)
        final_catalog.to_csv(csv_path, index=False)
        dbg(f"✅ Saved {csv_filename}")
        
        return True, len(final_catalog), msto_mag, msto_color
        
    except Exception as e:
        stage(f"❌ ERROR processing combination DM={apparent_distance_modulus:.2f} EBV={ebv:.2f}: {e}")
        if DEBUG:
            import traceback
            traceback.print_exc()
        return False, 0, np.nan, np.nan

def main():
    # Check command line arguments
    if len(sys.argv) != 3:
        stage("Usage: python scmd_grid_red_dist_chaboyer_g_err_3err.py <run_number> <age_gyr>")
        sys.exit(1)
    
    # Get run number and age from command line
    run_number = sys.argv[1]
    target_age_gyr = float(sys.argv[2])
    target_age_years = target_age_gyr * 1e9  # Convert to years
    
    # Construct file paths based on run number
    #isochrone_file = f"/dartfs-hpc/rc/lab/C/ChaboyerB/Isochrones_Only_July_10/{run_number}/isos/isochrones_example.txt"
    isochrone_file = f"/dartfs-hpc/scratch/f007058/Isochrones_W/{run_number}/isos/isochrones_example.txt"

    #/dartfs-hpc/scratch/f007058/Isochrones_W
    #output_dir = f"/dartfs-hpc/rc/lab/C/ChaboyerB/sCMD_July_11/{run_number}"
    output_dir = f"/dartfs-hpc/scratch/f007058/sCMD_W/{run_number}"
    #/dartfs-hpc/scratch/f007058/sCMD_W

    # Fixed file paths - updated to use regional scatter data
   # gaia_file = "/dartfs-hpc/rc/home/8/f007058/cleaned_ngc6791_no_blue_stragglers.csv"
    gaia_file = "/dartfs-hpc/rc/home/8/f007058/NGC6791_referee_isochrones_sCMDs/chaboyer_DM13.280_EBV0.200_age8.3Gyr_001.csv"
    completeness_file = "/dartfs-hpc/rc/home/8/f007058/g_mag_star_counts_and_completeness.csv"
    red_scatter_file = "/dartfs-hpc/rc/home/8/f007058/ngc6791_bin_statistics_red_region.csv"
    blue_scatter_file = "/dartfs-hpc/rc/home/8/f007058/ngc6791_bin_statistics_blue_region.csv"
    
    # Load regional scatter statistics
    red_stats, blue_stats = load_regional_scatter_data(red_scatter_file, blue_scatter_file)
    
    # Check if isochrone file exists
    if not os.path.exists(isochrone_file):
        dbg(f"❌ Isochrone file not found: {isochrone_file}")
        sys.exit(1)
    
    dbg(f"🚀 STARTING GRID SEARCH FOR RUN {run_number}")
    dbg(f"   Isochrone file: {isochrone_file}")
    dbg(f"   Output directory: {output_dir}")
    dbg(f"   Target age: {target_age_gyr:.1f} Gyr ({target_age_years:.0e} years)")
    
    # Get all available ages from the isochrone file
    all_ages_years = get_available_ages_from_isochrone(isochrone_file)
    
    # Find closest match to the specified age
    dbg(f"\n🎯 FINDING CLOSEST ISOCHRONE AGE:")
    dbg(f"   Target: {target_age_years:.0e} years ({target_age_gyr:.1f} Gyr)")
    
    closest_idx = np.argmin(np.abs(all_ages_years - target_age_years))
    closest_age = all_ages_years[closest_idx]
    
    if np.abs(closest_age - target_age_years) / target_age_years < 0.01:  # Within 1%
        selected_ages = [closest_age]
        dbg(f"✅ Found match: {closest_age:.0e} years (target: {target_age_years:.0e})")
    else:
        dbg(f"⚠️ No close match for {target_age_years:.0e} years (closest: {closest_age:.0e})")
        dbg(f"   Difference: {np.abs(closest_age - target_age_years)/target_age_years*100:.1f}%")
        # Use the closest available age anyway
        selected_ages = [closest_age]
        dbg(f"✅ Using closest available age: {closest_age:.0e} years")
    
    # Create specific parameter grid
    parameter_grid = create_parameter_grid_specific()
    
    # Calculate total combinations
    total_combinations = len(selected_ages) * len(parameter_grid)
    stage(f"[run {run_number} age {target_age_gyr:.1f}Gyr] {len(parameter_grid)} combos")

    # Track results
    all_results = []
    successful_runs = 0
    failed_runs = 0
    total_stars = 0

    # Process each combination
    for age_idx, age_years in enumerate(selected_ages):
        age_gyr = age_years / 1e9

        for param_idx, params in enumerate(parameter_grid):
            distance_modulus = params['apparent_distance_modulus']
            ebv = params['ebv']

            # Process single combination with regional scatter stats
            success, n_stars, msto_mag, msto_color = process_single_isochrone_grid_quota(
                age_years, distance_modulus, ebv,
                isochrone_file, gaia_file, red_stats, blue_stats,
                completeness_file, output_dir
            )

            # Record results with both age formats
            result = {
                'run_number': run_number,
                'target_age_gyr': target_age_gyr,
                'age_years': age_years,
                'age_gyr': age_gyr,
                'distance_modulus': distance_modulus,
                'ebv': ebv,
                'n_stars': n_stars,
                'msto_mag': msto_mag,
                'msto_color': msto_color,
                'status': 'SUCCESS' if success else 'FAILED',
                'filename_csv': generate_filename(age_years, distance_modulus, ebv, "csv")
            }
            all_results.append(result)

            if success:
                successful_runs += 1
                total_stars += n_stars
            else:
                failed_runs += 1
                stage(f"  FAIL combo {param_idx+1}/{len(parameter_grid)} DM={distance_modulus:.2f} EBV={ebv:.2f}")

        stage(f"  age {age_gyr:.1f}Gyr done: {successful_runs} ok / {failed_runs} fail")
    
    # Save comprehensive results summary
    results_df = pd.DataFrame(all_results)
    summary_path = os.path.join(output_dir, f"grid_search_summary_{run_number}_{target_age_gyr:.1f}Gyr.csv")
    results_df.to_csv(summary_path, index=False)
    
    # Final report
    stage(f"[run {run_number} age {target_age_gyr:.1f}Gyr] DONE: {successful_runs}/{total_combinations} ok, {total_stars:,} stars")
    dbg(f"\n{'='*80}")
    dbg(f"🎉 GRID SEARCH COMPLETE FOR RUN {run_number}!")
    dbg(f"{'='*80}")
    dbg(f"✅ Target age: {target_age_gyr:.1f} Gyr")
    dbg(f"✅ Actual age used: {age_gyr:.1f} Gyr ({age_years:.0e} years)")
    dbg(f"✅ Successful runs: {successful_runs}/{total_combinations}")
    dbg(f"❌ Failed runs: {failed_runs}/{total_combinations}")
    dbg(f"📊 Total stars generated: {total_stars:,}")
    dbg(f"📁 Results saved in: {output_dir}/")
    dbg(f"📋 Summary file: grid_search_summary_{run_number}_{target_age_gyr:.1f}Gyr.csv")
    
    if successful_runs > 0:
        successful_results = results_df[results_df['status'] == 'SUCCESS']
        dbg(f"\n✅ Successfully processed combinations:")
        dbg(f"   Distance modulus range: {successful_results['distance_modulus'].min():.2f} - {successful_results['distance_modulus'].max():.2f}")
        dbg(f"   E(B-V) range: {successful_results['ebv'].min():.3f} - {successful_results['ebv'].max():.3f}")
        dbg(f"   Average stars per run: {successful_results['n_stars'].mean():.0f}")
        dbg(f"   Star count range: {successful_results['n_stars'].min()} - {successful_results['n_stars'].max()}")
    
    dbg(f"{'='*80}")

def get_available_ages_from_isochrone(iso_path):
    """
    Extract all available ages from the isochrone file.
    
    Parameters:
    -----------
    iso_path : str
        Path to the isochrone file
        
    Returns:
    --------
    ages_years : np.ndarray
        Array of available ages in years
    """
    # Load isochrone file (cached) and get unique ages
    iso_df = _load_iso_raw(iso_path)
    unique_ages = np.unique(iso_df["isochrone_age_yr"])
    
    dbg(f"✅ Found {len(unique_ages)} unique ages in isochrone file:")
    for age in unique_ages:
        age_gyr = age / 1e9
        dbg(f"   {age:.0e} years ({age_gyr:.1f} Gyr)")
    
    return unique_ages

def create_parameter_grid_specific():
    """
    Create a specific parameter grid for apparent distance modulus and extinction.
    
    Returns:
    --------
    parameter_grid : list
        List of parameter combinations
    """
    # Define parameter ranges - these are APPARENT distance moduli from literature
    apparent_distance_moduli = np.arange(13.23, 13.51, 0.01)  # Apparent DM from literature
    ebv_values = np.arange(0.10, 0.21, 0.01)                 # E(B-V) values
    
    dbg(f"🔧 Creating parameter grid:")
    dbg(f"   Apparent distance modulus: {apparent_distance_moduli.min():.2f} to {apparent_distance_moduli.max():.2f}")
    dbg(f"   E(B-V): {ebv_values.min():.2f} to {ebv_values.max():.2f}")
    
    # Create all combinations
    parameter_grid = []
    combination_id = 1
    
    for dm in apparent_distance_moduli:
        for ebv in ebv_values:
            parameter_grid.append({
                'combination_id': combination_id,
                'apparent_distance_modulus': dm,  # Changed name to be clear
                'ebv': ebv
            })
            combination_id += 1
    
    dbg(f"   Total combinations: {len(parameter_grid)}")
    
    return parameter_grid

def generate_filename(age_years, distance_modulus, ebv, extension):
    """
    Generate standardized filenames for output files.
    
    Parameters:
    -----------
    age_years : float
        Age in years
    distance_modulus : float
        Distance modulus
    ebv : float
        E(B-V) extinction
    extension : str
        File extension (csv, png, etc.)
        
    Returns:
    --------
    filename : str
        Standardized filename
    """
    age_gyr = age_years / 1e9
    filename = f"NGC6791_age{age_gyr:.1f}Gyr_DM{distance_modulus:.3f}_EBV{ebv:.3f}.{extension}"
    return filename

def quota_sample_stars(iso_df_out, gaia_file, completeness_file, imf_slope=-1.02):
    """
    Generate synthetic stars using completeness-corrected sampling.
    
    Steps:
    1. Count stars per bin from observational data (G mag 16.25-18.0)
    2. Correct for completeness by dividing by completeness fraction
    3. Sum to get total stars to generate
    4. Generate EXACTLY that many stars using IMF sampling in mass range for G 16.25-18.0
    5. Apply completeness by MULTIPLYING generated count by completeness per bin
    """
    # Load observational data
    gaia_data = _GAIA_CACHE.get(gaia_file)
    if gaia_data is None:
        gaia_data = pd.read_csv(gaia_file)
        _GAIA_CACHE[gaia_file] = gaia_data
    completeness_data = _COMPLETENESS_CACHE.get(completeness_file)
    if completeness_data is None:
        completeness_data = pd.read_csv(completeness_file)
        _COMPLETENESS_CACHE[completeness_file] = completeness_data
    
    dbg(f"✅ Loaded observational data: {len(gaia_data)} stars")
    dbg(f"✅ Loaded completeness data: {len(completeness_data)} bins")
    
    # Define bin edges from completeness file (0.5 mag bins)
    bin_centers = completeness_data['bin_center'].values
    completeness_fractions = completeness_data['completeness'].values
    
    # Calculate bin edges
    bin_edges = []
    for center in bin_centers:
        bin_edges.append(center - 0.25)
    bin_edges.append(bin_centers[-1] + 0.25)
    bin_edges = np.array(bin_edges)
    
    dbg(f"\n🔧 Bin structure:")
    dbg(f"   Bin edges: {bin_edges[0]:.2f} to {bin_edges[-1]:.2f}")
    dbg(f"   Number of bins: {len(bin_centers)}")
    
    # ⭐ EXACT RANGE: 16.25-18.0 ⭐
    G_MIN = 16.25
    G_MAX = 18.0
    
    # Filter observational data to G range
    g_obs = gaia_data['Gmag'].values
    obs_mask = (g_obs >= G_MIN) & (g_obs <= G_MAX)
    g_obs_filtered = g_obs[obs_mask]
    
    dbg(f"\n🔧 Observational data (G {G_MIN}-{G_MAX}): {len(g_obs_filtered)} stars")
    
    # Count stars per bin from observational data
    observed_counts = []
    completeness_corrected_counts = []
    valid_bin_centers = []
    valid_completeness = []
    
    for i, center in enumerate(bin_centers):
        # Only process bins in our G range
        if center < G_MIN or center > G_MAX:
            continue
            
        left = bin_edges[i]
        right = bin_edges[i + 1]
        
        # Count observed stars in this bin
        bin_mask = (g_obs_filtered >= left) & (g_obs_filtered < right)
        n_observed = np.sum(bin_mask)
        
        # Correct for completeness
        if completeness_fractions[i] > 0:
            n_corrected = n_observed / completeness_fractions[i]
        else:
            n_corrected = 0
        
        observed_counts.append(n_observed)
        completeness_corrected_counts.append(n_corrected)
        valid_bin_centers.append(center)
        valid_completeness.append(completeness_fractions[i])
        
        dbg(f"   Bin {center:.2f}: {n_observed} observed → {n_corrected:.1f} corrected (completeness={completeness_fractions[i]:.3f})")
    
    # Total stars to generate
    total_stars_needed = int(np.sum(completeness_corrected_counts))
    
    dbg(f"\n✅ Total stars to generate: {total_stars_needed}")
    
    # ⭐ FIND EXACT MASS RANGE - NO BUFFER ⭐
    mass_iso = iso_df_out["initial_mass"].values
    G_iso = iso_df_out["Gaia_G_EDR3"].values
    
    # Sort by mass
    sort_idx = np.argsort(mass_iso)
    m_sorted = mass_iso[sort_idx]
    G_sorted = G_iso[sort_idx]
    
    # Find exact mass range for G 16.25-18.0
    g_mask = (G_sorted >= G_MIN) & (G_sorted <= G_MAX)
    
    if np.sum(g_mask) == 0:
        dbg(f"❌ No isochrone points in G range {G_MIN}-{G_MAX}!")
        return pd.DataFrame()
    
    m_min_for_g_range = m_sorted[g_mask].min()
    m_max_for_g_range = m_sorted[g_mask].max()
    
    dbg(f"\n🔧 Mass range for G {G_MIN}-{G_MAX}:")
    dbg(f"   Mass: {m_min_for_g_range:.3f} to {m_max_for_g_range:.3f} M☉")
    
    # Create interpolators
    BP_iso = iso_df_out["Gaia_BP_EDR3"].values[sort_idx]
    RP_iso = iso_df_out["Gaia_RP_EDR3"].values[sort_idx]
    
    G_interp = interp1d(m_sorted, G_sorted, kind="linear", fill_value="extrapolate")
    BP_interp = interp1d(m_sorted, BP_iso, kind="linear", fill_value="extrapolate")
    RP_interp = interp1d(m_sorted, RP_iso, kind="linear", fill_value="extrapolate")
    
    # IMF sampling function - ONLY in the restricted mass range
    def sample_imf(n, m_min, m_max, slope):
        r = np.random.random(n)
        exponent = slope + 1.0
        return (r * (m_max ** exponent - m_min ** exponent) + m_min ** exponent) ** (1/exponent)
    
    # ⭐ SAMPLE EXACTLY THE NUMBER NEEDED - NO OVERSAMPLING ⭐
    dbg(f"\n🔧 Sampling {total_stars_needed} stars in mass range {m_min_for_g_range:.3f}-{m_max_for_g_range:.3f} M☉...")
    
    primary_masses = sample_imf(total_stars_needed, m_min_for_g_range, m_max_for_g_range, imf_slope)
    
    # Interpolate magnitudes
    sim_G = G_interp(primary_masses)
    sim_BP = BP_interp(primary_masses)
    sim_RP = RP_interp(primary_masses)
    
    dbg(f"✅ Generated {len(sim_G)} stars")
    dbg(f"   G range: {sim_G.min():.2f} to {sim_G.max():.2f}")
    
    # ⭐ NOW APPLY COMPLETENESS BY MULTIPLYING GENERATED COUNT ⭐
    dbg(f"\n🔧 Applying completeness by multiplying generated count per bin...")
    
    keep_mask = np.ones(len(sim_G), dtype=bool)
    
    for bin_idx, center in enumerate(valid_bin_centers):
        # Calculate bin edges for this valid bin
        left = center - 0.25
        right = center + 0.25
        
        # Find stars in this bin
        bin_mask = (sim_G >= left) & (sim_G < right)
        bin_indices = np.where(bin_mask)[0]
        n_in_bin = len(bin_indices)
        
        if n_in_bin > 0:
            # ⭐ MULTIPLY BY COMPLETENESS to get final count ⭐
            completeness = valid_completeness[bin_idx]
            n_to_keep = int(round(n_in_bin * completeness))
            n_to_remove = n_in_bin - n_to_keep
            
            if n_to_remove > 0:
                # Randomly select stars to remove
                remove_indices = np.random.choice(bin_indices, size=n_to_remove, replace=False)
                keep_mask[remove_indices] = False
            
            dbg(f"   Bin {center:.2f}: {n_in_bin} generated × {completeness:.3f} = {n_to_keep} kept (removing {n_to_remove})")
    
    # Apply the keep mask
    sim_G = sim_G[keep_mask]
    sim_BP = sim_BP[keep_mask]
    sim_RP = sim_RP[keep_mask]
    primary_masses = primary_masses[keep_mask]
    
    # Create final catalog
    synthetic_catalog = pd.DataFrame({
        'mass': primary_masses,
        'G': sim_G,
        'BP': sim_BP,
        'RP': sim_RP
    })
    
    dbg(f"\n✅ Final catalog: {len(synthetic_catalog)} stars after completeness correction")
    if len(synthetic_catalog) > 0:
        dbg(f"   G magnitude range: {synthetic_catalog['G'].min():.2f} to {synthetic_catalog['G'].max():.2f}")
    
    return synthetic_catalog

def calculate_perp_distance_to_ridge_normalized(star_colors, star_mags, ridge_colors, ridge_mags, 
                                                 color_scale, mag_scale, n_interp_points=1000):
    """
    Calculate perpendicular distance with PARAMETRIC interpolation to handle vertical sections.
    """
    from scipy.interpolate import splprep, splev
    
    # Normalize coordinates
    star_colors_norm = star_colors / color_scale
    star_mags_norm = star_mags / mag_scale
    ridge_colors_norm = ridge_colors / color_scale
    ridge_mags_norm = ridge_mags / mag_scale
    
    # ⭐ USE PARAMETRIC SPLINE (handles vertical/horizontal sections) ⭐
    # This creates a parameter t that goes from 0 to 1 along the curve
    # instead of trying to express color as a function of magnitude
    
    # Sort by magnitude first (but parametric spline doesn't require monotonic relationship)
    sort_idx = np.argsort(ridge_mags_norm)
    ridge_colors_sorted = ridge_colors_norm[sort_idx]
    ridge_mags_sorted = ridge_mags_norm[sort_idx]
    
    # Create parametric spline: both color and mag are functions of parameter t
    tck, u = splprep([ridge_colors_sorted, ridge_mags_sorted], s=0, k=3)
    
    # Evaluate spline at many points
    u_fine = np.linspace(0, 1, n_interp_points)
    ridge_colors_dense, ridge_mags_dense = splev(u_fine, tck)
    
    dbg(f"   Interpolated ridge from {len(ridge_colors_norm)} to {len(ridge_colors_dense)} points using parametric spline")
    
    # Create (N_ridge x 2) array of dense ridge points
    ridge_points = np.column_stack([ridge_colors_dense, ridge_mags_dense])
    
    # Create (N_stars x 2) array of star positions
    star_points = np.column_stack([star_colors_norm, star_mags_norm])
    
    # For each star, find minimum distance to any ridge segment
    distances = np.zeros(len(star_colors))
    
    for i in range(len(star_points)):
        star = star_points[i]
        
        # Distance to all ridge points
        point_dists = np.sqrt(np.sum((ridge_points - star)**2, axis=1))
        
        # Find closest ridge point
        min_idx = np.argmin(point_dists)
        
        # Check segments around closest point (check a few neighbors)
        check_start = max(0, min_idx - 2)
        check_end = min(len(ridge_points) - 1, min_idx + 3)
        check_indices = range(check_start, check_end)
        
        min_seg_dist = np.inf
        for j in check_indices:
            if j >= len(ridge_points) - 1:
                break
                
            p1 = ridge_points[j]
            p2 = ridge_points[j + 1]
            
            # Project star onto segment
            v = p2 - p1
            w = star - p1
            
            seg_len_sq = np.dot(v, v)
            if seg_len_sq < 1e-10:
                dist = np.linalg.norm(w)
            else:
                t = np.clip(np.dot(w, v) / seg_len_sq, 0, 1)
                projection = p1 + t * v
                dist = np.linalg.norm(star - projection)
            
            min_seg_dist = min(min_seg_dist, dist)
        
        distances[i] = min_seg_dist
    
    return distances

def apply_scatter_to_region(region_G, region_BP, region_RP, region_stats, slope_interpolator, 
                            color_scale, mag_scale, region_name, 
                            max_perp_distance_normalized=0.05):
    """
    Apply scatter to stars in a specific region using that region's statistics.
    ⭐ CLIPS scatter to max normalized perpendicular distance ⭐
    """
    n_stars = len(region_G)
    if n_stars == 0:
        return {'G': np.array([]), 'BP': np.array([]), 'RP': np.array([])}
    
    dbg(f"   Processing {region_name} region: {n_stars} stars")
    
    # Get bin information from region_stats
    bin_lefts = region_stats['bin_left'].values
    bin_rights = region_stats['bin_right'].values
    perp_std_vals = region_stats['perp_std'].values  # In PHYSICAL units
    
    dbg(f"     Regional bins: {len(region_stats)} bins")
    dbg(f"     G range: {bin_lefts.min():.2f} to {bin_rights.max():.2f}")
    dbg(f"     Perp_std range: {perp_std_vals.min():.4f} to {perp_std_vals.max():.4f} (physical units)")
    
    # Vectorized bin assignment for all stars based on G magnitude
    bin_assignments = np.full(n_stars, -1, dtype=int)
    
    for bin_idx in range(len(region_stats)):
        g_left = bin_lefts[bin_idx]
        g_right = bin_rights[bin_idx]
        mask = (region_G >= g_left) & (region_G < g_right)
        bin_assignments[mask] = bin_idx
    
    # Only process stars with valid bins
    valid_mask = bin_assignments >= 0
    n_valid = np.sum(valid_mask)
    
    if n_valid == 0:
        return {'G': np.array([]), 'BP': np.array([]), 'RP': np.array([])}
    
    # Get valid star data
    valid_g = region_G[valid_mask]
    valid_bp = region_BP[valid_mask]
    valid_rp = region_RP[valid_mask]
    valid_bins = bin_assignments[valid_mask]
    
    # Initialize output arrays for valid stars
    sim_G_obs = valid_g.copy()
    sim_BP_obs = valid_bp.copy()
    sim_RP_obs = valid_rp.copy()
    
    # Get isochrone slopes for all valid stars
    iso_slopes = slope_interpolator(valid_g)
    
    # Calculate perpendicular slopes
    perp_slopes = np.where(np.abs(iso_slopes) > 1e-6, 
                          -1.0 / iso_slopes, 
                          np.sign(iso_slopes) * 1e6)
    
    # Calculate unit perpendicular vectors
    perp_norms = np.sqrt(1 + perp_slopes**2)
    unit_perp_g = 1.0 / perp_norms
    unit_perp_color = perp_slopes / perp_norms;
    
    # Get sigma values for each star
    sigma_scatter_vals = perp_std_vals[valid_bins]
    
    # Sample error in PHYSICAL space
    error_magnitudes_physical = np.random.normal(0, sigma_scatter_vals)
    
    # ⭐ CLIP ERROR MAGNITUDES TO MAX NORMALIZED DISTANCE ⭐
    # Calculate what this error would be in normalized space
    g_errors_unnormalized = error_magnitudes_physical * unit_perp_g
    color_errors_unnormalized = error_magnitudes_physical * unit_perp_color
    
    # Normalize
    g_errors_norm = g_errors_unnormalized / mag_scale
    color_errors_norm = color_errors_unnormalized / color_scale
    
    # Calculate normalized distance
    distances_norm = np.sqrt(g_errors_norm**2 + color_errors_norm**2)
    
    # Clip to max distance
    scale_factors = np.where(distances_norm > max_perp_distance_normalized,
                            max_perp_distance_normalized / distances_norm,
                            1.0)
    
    # Apply scale factors to errors
    g_errors = g_errors_unnormalized * scale_factors
    color_errors = color_errors_unnormalized * scale_factors
    
    # Count clipped stars
    n_clipped = np.sum(scale_factors < 1.0)
    
    # Apply errors
    sim_G_obs += g_errors
    
    bp_errors = color_errors * 0.5
    rp_errors = -color_errors * 0.5;
    
    sim_BP_obs += bp_errors
    sim_RP_obs += rp_errors
    
    dbg(f"     ✅ Applied scatter to {n_valid} stars in {region_name} region")
    if n_clipped > 0:
        dbg(f"     ⚠️ Clipped {n_clipped} stars ({100*n_clipped/n_valid:.1f}%) to max distance")
    dbg(f"     Final G range: {sim_G_obs.min():.2f} to {sim_G_obs.max():.2f}")
    
    return {
        'G': sim_G_obs,
        'BP': sim_BP_obs,
        'RP': sim_RP_obs
    }


def add_regional_g_magnitude_scatter(sim_G, sim_BP, sim_RP, iso_df_out, red_stats, blue_stats,
                                    max_perp_distance_normalized=0.05):
    """
    Add perpendicular scatter with clipping - NO additional filtering needed!
    """
    dbg("🔧 Adding regional G magnitude-based scatter perpendicular to isochrone...")
    
    # Convert to numpy arrays
    sim_G = np.array(sim_G, dtype=float)
    sim_BP = np.array(sim_BP, dtype=float) 
    sim_RP = np.array(sim_RP, dtype=float)
    sim_color = sim_BP - sim_RP
    
    # ⭐ USE FIXED DATA BOUNDARIES FROM NGC 6791 REFERENCE ⭐
    bp_rp_min = 0.9020960000000002
    bp_rp_max = 1.7554049999999997
    g_min = 15.588097
    g_max = 18.999708
    
    color_scale = bp_rp_max - bp_rp_min
    mag_scale = g_max - g_min
    
    dbg(f"   Using FIXED NGC 6791 data boundaries for normalization:")
    dbg(f"   BP-RP range: {bp_rp_min:.6f} to {bp_rp_max:.6f} (scale={color_scale:.6f})")
    dbg(f"   G range: {g_min:.6f} to {g_max:.6f} (scale={mag_scale:.6f})")
    dbg(f"   Max perpendicular distance: {max_perp_distance_normalized} (normalized units)")
    
    # Get isochrone for filtering
    color_iso = iso_df_out["Gaia_BP_EDR3"] - iso_df_out["Gaia_RP_EDR3"]
    g_iso = iso_df_out["Gaia_G_EDR3"]
    
    iso_mask = (g_iso >= g_min - 0.5) & (g_iso <= g_max + 0.5)
    color_iso_filtered = color_iso[iso_mask].values
    g_iso_filtered = g_iso[iso_mask].values
    
    dbg(f"   Filtered isochrone: {len(color_iso_filtered)} points")
    
    # Calculate isochrone slopes
    slope_interpolator = calculate_isochrone_slope_g_vs_color(iso_df_out)
    
    # Classify stars into regions
    red_region_mask = (sim_color >= 1.15) & (sim_G < 17.5)
    blue_region_mask = ~red_region_mask
    
    dbg(f"   Red region stars: {np.sum(red_region_mask)}")
    dbg(f"   Blue region stars: {np.sum(blue_region_mask)}")
    
    # Apply scatter to each region (WITH CLIPPING)
    all_G = []
    all_BP = []
    all_RP = []
    
    # Process red region
    if np.sum(red_region_mask) > 0:
        red_applied = apply_scatter_to_region(
            sim_G[red_region_mask], sim_BP[red_region_mask], sim_RP[red_region_mask],
            red_stats, slope_interpolator, color_scale, mag_scale, "RED",
            max_perp_distance_normalized
        )
        if len(red_applied['G']) > 0:
            all_G.extend(red_applied['G'])
            all_BP.extend(red_applied['BP'])
            all_RP.extend(red_applied['RP'])
    
    # Process blue region
    if np.sum(blue_region_mask) > 0:
        blue_applied = apply_scatter_to_region(
            sim_G[blue_region_mask], sim_BP[blue_region_mask], sim_RP[blue_region_mask],
            blue_stats, slope_interpolator, color_scale, mag_scale, "BLUE",
            max_perp_distance_normalized
        )
        if len(blue_applied['G']) > 0:
            all_G.extend(blue_applied['G'])
            all_BP.extend(blue_applied['BP'])
            all_RP.extend(blue_applied['RP'])
    
    # Convert to arrays
    sim_G_obs = np.array(all_G)
    sim_BP_obs = np.array(all_BP)
    sim_RP_obs = np.array(all_RP)
    sim_color_obs = sim_BP_obs - sim_RP_obs
    
    # ⭐ NO FILTERING - clipping already done! ⭐
    dbg(f"✅ Applied regional scatter with clipping")
    dbg(f"   Total stars: {len(sim_G_obs)}")
    if len(sim_G_obs) > 0:
        dbg(f"   Final G range: {sim_G_obs.min():.2f} to {sim_G_obs.max():.2f}")
        dbg(f"   Final BP-RP range: {sim_color_obs.min():.2f} to {sim_color_obs.max():.2f}")
    
    return sim_G_obs, sim_BP_obs, sim_RP_obs, sim_color_obs


if __name__ == "__main__":
    main()