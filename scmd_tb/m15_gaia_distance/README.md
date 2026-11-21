# M15: do the Gaia distance or Ying's own (mu, E) close the age gap? (2026-11-11)

Purpose: test whether an external distance (Gaia), or Ying et al. (2025)'s own best-fit distance and reddening,
moves the production M15 age (14.2 Gyr) toward Ying's 13.39 +/- 0.39 Gyr. No new grid was run. The production
setup and code are unchanged.

Conventions: the grid's mu is the apparent F606W modulus, and E is the F606W-F814W shift (reddening "evi"). A
true modulus therefore means mu = mu0 + 2.82 E (ACS: A_F606W 2.88 E(B-V), E(F606W-F814W) 1.02 E(B-V)).
Baumgardt & Vasiliev (2021): 10.71 +/- 0.10 kpc, i.e. mu0 = 15.149 +/- 0.020. Ying's Table 7 (M15, 2D KS): DM
15.49 +/- 0.03, E 0.08 +/- 0.01, [Fe/H] -2.30, so mu0 = 15.264.

## Results

1. **Every free fit implies mu0 of about 15.30, at every age** (`logs/ying_point_check.out`). The best fit at
   each age from 12.0 to 14.6 Gyr gives mu0 = 15.29-15.37 in production, and the same with the map off. Ying's
   own fit gives 15.26. All of these are 0.11-0.22 mag beyond Gaia. Age trades against E along a line of
   constant mu0, so a distance prior cannot choose between ages. It disagrees with all of them, Ying's
   included. That points to an absolute offset shared by both pipelines (photometric zero points, bolometric
   corrections, or the model main sequence), not to something that separates us from Ying.
2. **At Ying's own distance and reddening, our pipeline still gives 14.2 Gyr:**

| cells restricted to | production (map on) | map off |
|---|---|---|
| mu 15.46-15.52, E 0.07-0.09 (Ying +/- 1 sigma) | 14.2 Gyr, KS 0.0139 (13.4 Gyr: 0.0230) | 14.1 Gyr, KS 0.0192 (13.4 Gyr: 0.0248) |
| mu 15.43-15.55, E 0.06-0.10 (+/- 2 sigma) | 14.2 Gyr, KS 0.0139 (13.4 Gyr: 0.0153) | 14.0 Gyr, KS 0.0190 (13.4 Gyr: 0.0246) |

   Production's global best is 0.0137, so a 14.2 Gyr fit sits essentially at Ying's point. The age gap is not
   in the distance or the reddening.
3. **The scan along the Gaia line is inconclusive** (`gaia_line_scan.py`, `logs/gaia_line_scan.out`,
   `gaia_line_scan.csv`). It covered runs 10357, 10392 and 10245, ages 12.6-14.6 Gyr, mu0 15.149 and 15.189,
   and E 0.06-0.16. Every value came out at KS 0.036-0.19, far above the free fits. But the E step (0.01) is
   coarser than the KS valley (0.001 in E already moves the KS from 0.0137 to 0.0196), so the minima along the
   line were missed. Do not use its per-age numbers. Among the existing grid cells within 0.04 of Gaia's mu0,
   the lowest KS is at 14.3-14.6 Gyr (0.016-0.018). Those cells were not optimised for the constraint, so this
   is suggestive only.

## Files

- `gaia_line_scan.py` and `logs/gaia_line_scan.out`, `gaia_line_scan.csv`: the scan in result 3.
- `ying_point_check.py` and `logs/ying_point_check.out`: results 1 and 2. `ying_point_check.py` filters
  `../m15_final_refined/combined_2d_ks_results.csv` and `../m15_yingbox_nozp_refined/combined_2d_ks_results.csv`
  by mu and E and takes the minimum KS per age, with ages capped at each run's age limit.
