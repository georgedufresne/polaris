# M15: calibration-star weights (2026-11-11)

Purpose: apply Ying et al.'s calibration-star weighting to our 500 M15 Monte Carlo runs and see whether it moves
the age. The production fits are unchanged. Only the posterior weights change.

## Recipe (from the papers)

- **Stars:** Ying et al. 2023 (M92), section 2.1: HIP 46120, M_F606W 5.7867 +/- 0.0026, F606W-F814W
  0.566 +/- 0.002. HIP 106924, 6.0406 +/- 0.0037, 0.601 +/- 0.005. Gaia EDR3 parallaxes with half the
  Lindegren zero-point correction, zero reddening. Ying et al. 2025 (Table 2) lists the same two stars for M15
  ([Fe/H] -2.22 and -2.23) and scores them against the cluster's own Monte Carlo isochrones.
- **Weight:** 2025 Eq. 4-5 and 2023 Eq. 1-2. For each age, chi2 = sum over the two stars of
  ((V_iso - V_star)/sigma_V)^2 + ((VI_iso - VI_star)/sigma_VI)^2. Each run takes its lowest chi2 over ages. Its
  weight is w_cal = 1 - CDF_chi2,2dof = exp(-chi2_min/2), and the final weight is w_cal x w_CMD.
- **Our implementation:** the squared normalised distance from each star to the isochrone curve (linear between
  EEPs), over ages 8-16 Gyr capped at each run's age limit. The isochrones and MIST ACS/WFC bolometric
  corrections are the same as in production. The papers do not say which isochrone point is compared. Along the
  main sequence, the normalised distance and the colour difference at the star's magnitude differ by about 3%.

## Results (`calib_weights.py`, `logs/calib_weights.out`, `calib_chi2.csv`)

- **Our isochrones pass the test easily.** chi2_min has a median of 0.7. 59% of runs have chi2 < 2 (the M92
  paper: about 20% of its sets), and 67% have p > 0.05 (M92: about 35%). The median colour offset at the stars'
  magnitudes is +0.0002 mag. The offset tracks mixing length (Spearman -0.75 with CMIXLA) and almost nothing
  else.
- **The weights barely move the age** (mean, median, 16-84% in Gyr):

| posterior | without the weights | with the weights |
|---|---|---|
| production strict | 14.22, 14.20, 14.10-14.40 | 14.20, 14.20, 14.10-14.30 |
| production shifted | 14.25, 14.30, 14.10-14.40 | 14.26, 14.30, 14.10-14.50 |
| production scaled | 14.10, 14.20, 13.80-14.40 | 14.13, 14.20, 13.90-14.50 |
| no map, Ying's box, shifted | 13.89, 14.20, 13.10-14.40 | 13.80, 14.10, 12.90-14.40 |
| no map, Ying's box, scaled | 13.15, 13.20, 12.20-14.20 | 13.08, 13.00, 12.20-14.10 |

- **The best-fitting runs at every age already pass the stars** (p > 0.05). Best production fit: 14.2 Gyr,
  KS 0.0139 among passing runs (0.0137 over all runs).
- **A side result:** the models match stars of known parallax on the lower main sequence. Yet every M15 fit
  implies mu0 of about 15.30 against Gaia's 15.15 (`../m15_gaia_distance/`), and Ying's own fit implies 15.26.
  The distance tension is therefore probably on the M15 side (reddening, the photometric zero points between the
  two HST datasets, or the cluster distance), not in the isochrones' lower main sequence. This is an inference.

Conclusion: calibration-star weighting does not explain the gap. It shifts the production age by 0.03 Gyr or
less.

## Caveats

- The bolometric corrections are MIST at [Fe/H] -2.25 for every run, as in production. Ying's bolometric
  corrections are not stated (the 2025 paper mentions PHOENIX models at the nearest [alpha/Fe]). The test
  validates our isochrones together with our bolometric corrections, which is the combination we fit with.
- The stars are main-sequence dwarfs. They test the lower main sequence, mostly through mixing length, not the
  turnoff or subgiant branch that set the age.
