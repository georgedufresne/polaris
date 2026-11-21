# M15 recreation with Ying's published isochrones (2026-11-11)

Purpose: recreate Ying et al. (2025)'s M15 fit as far as his logic is known, with his own isochrones, keeping our
sCMD construction and our 2D KS. It is a one-change test against `../m15_yingbox_nozp` (same setup, our DSEP4
isochrones).

| | this run | `../m15_yingbox_nozp` |
|---|---|---|
| isochrones | Ying's, Zenodo 10.5281/zenodo.13998637 (his absolute F606W and F606W-F814W) | our DSEP4 rebuild + MIST ACS/WFC BCs |
| ages | 8.0-16.0 Gyr in 0.2 (his grid) | 8.0-15.0 in 0.1 |
| zero-point map | off (as Ying) | off |
| distance and reddening | mu 15.30-15.60, E 0.08-0.15 applied directly (as Ying) | same |
| sCMDs, KS, refinement, bootstrap, posteriors | ours: ridge cut, nearest-AST errors, PDMF -1.00, binaries 0.017, 95th-percentile KS | same |

Same 500 runs. His input parameters equal ours for every run except core overshoot, which ours set to 0 (see
`../m15_ying_isochrones/`).

## Files

- **`convert_ying_isochrones.py`** (`logs/convert_ying_isochrones.out`) writes his HDF5 isochrones for our 500 runs
  as `/dartfs-hpc/scratch/f007058/M15/iso_ying2025/<run>/ying_cmd.<run>`, in the DSEP block format the pipeline
  already reads (`iso_format` "dsep"), and writes `age_limits_ying.json` there. All 500 read back exactly, and
  every run is valid to 16 Gyr.
- **`setup.json`** is `../m15_yingbox_nozp/setup.json` with `iso`, `iso_format` and `ages` changed. The other
  inputs:
  - `refine/setup_override.json`: E 0.080-0.150 in 0.001, as there;
  - `code/`: a copy of `../m15_yingbox_nozp/code`, unchanged;
  - `runs.txt`: the same 500 runs.
- **`m15yi_statistics.py`** is `../m15f_statistics.py` with the histogram bins widened to 16 Gyr (0.2 steps) and
  E 0.155. The posterior summaries are computed from these histograms, so rows outside them would be dropped.
- **`validate_yingiso.py`** and `validate.sbatch` (`logs/validate_11026029.out`).
- **`run.sh`** (driver, `driver.pid`, `run.log`): grid (`m15yi_scmd`), collect, KS refinement (`m15yi_refine`),
  merge, then `../m15_yingiso_nozp_refined/`: bootstrap (`boot_m15yi`), statistics, plots, and
  `compare_ks_vs_age.py` (`analysis/ks_vs_age_our_vs_ying_isochrones.png`).

## Validation (2026-11-11)

1. This code copy with `../m15_yingbox_nozp_refined/setup.json` reproduces that run's best cell (10245, 14.0 Gyr,
   mu 15.54, E 0.085): KS 0.0189500 vs 0.0189500 in its file.
2. With Ying's isochrones, the same cell gives KS 0.0893 (32,089 fit stars, map off). His 14.0 Gyr isochrone for
   this run is not ours, so the minimum moves.
3. At Ying's published best fit (DM 15.49, E 0.08), KS with his isochrones / ours:

| run | 13.4 Gyr | 14.2 Gyr |
|---|---|---|
| 10245 | 0.0237 / 0.1170 | 0.0658 / 0.0252 |
| 10357 | 0.1967 / 0.1884 | 0.2843 / 0.1240 |
| 10392 | 0.5360 / 0.1531 | 0.6303 / 0.2610 |

   His isochrones fit reasonably near his published point for some runs (10245 at 13.4 Gyr). The very large values
   for 10392 come from his isochrone being about 0.05 redder at the turnoff (`../m15_ying_isochrones/`), which
   E 0.08 cannot absorb.

## Results (finished 2026-11-12 00:33)

The run went grid 2026-11-11 14:12-22:49 (2,500 units) -> refinement to 00:24 (0 errors) -> bootstrap at the
best cell to 00:33 -> statistics. Outputs are in `../m15_yingiso_nozp_refined/` (posteriors, `analysis/`) and
`analysis/ks_vs_age_our_vs_ying_isochrones.png`.

- **Fit quality:**
  - Best KS 0.01341 (run 10089, 12.6 Gyr, mu 15.50, E 0.086), 2.2 sigma above its bootstrap null
    (0.00767 +/- 0.00262; 3 sigma line 0.01552). That is inside the null, and as good as production with the
    map (0.0137).
  - Next best: 10213 at 12.0 Gyr (0.01355) and 10027 at 11.6 Gyr (0.01364).
  - With our isochrones and the same setup the best is 0.0190, 4.3 sigma above the null.
- **Best KS by age:**

| age (Gyr) | 11.2 | 11.6 | 12.0 | 12.4 | 12.8 | 13.2 | 13.6 | 14.0 | 14.4 | 14.8 |
|---|---|---|---|---|---|---|---|---|---|---|
| Ying's isochrones, no map, his box (this run) | 0.0142 | 0.0136 | 0.0135 | 0.0138 | 0.0140 | 0.0138 | 0.0156 | 0.0175 | 0.0179 | 0.0198 |
| our isochrones, same setup | 0.0250 | 0.0234 | 0.0252 | 0.0226 | 0.0221 | 0.0214 | 0.0233 | 0.0190 | 0.0197 | 0.0369 |
| production (our isochrones, map, E 0-0.12) | 0.0277 | 0.0244 | 0.0262 | 0.0201 | 0.0184 | 0.0161 | 0.0150 | 0.0142 | 0.0147 | 0.0187 |

  With his isochrones the curve has a flat floor from 11.2 to 13.2 Gyr. It rises from 13.6 Gyr (outside 3 sigma
  from about 13.6 up), so ages of 13.6 Gyr and above are disfavoured.
- **Posteriors (age in Gyr):**

| | strict | shifted | scaled |
|---|---|---|---|
| KS weights | 12.21 (median 12.10, 16-84% 11.60-12.66) | 12.31 (11.55-12.88) | 12.39 (11.50-13.07) |
| times calibration weights (his isochrones; `calib_reweight.py`, `logs/calib_reweight.out`) | 12.15 (11.60-12.60) | 12.11 (11.60-12.80) | 12.09 (11.40-12.80) |

  Other strict medians: DM 15.53, E 0.087 (16-84% 0.080-0.095), [Fe/H] -2.28, mixing length 1.92.
- **Ying's M92 (2023) weighting rule** (`m92_rule_posterior.py`, `logs/m92_rule_posterior.out`, 2026-11-14): one best
  fit per run, keep runs within 3 sigma of the null mean, fit probability 1 below the mean and from the null's tail
  above it, times calibration weights. 119 of 500 runs pass; age 12.21 (12.16 with calibration weights), the same as
  our strict posterior. Both readings of his unstated formula give the same age. Applied to `../m15_yingiso_err99`,
  no run passes. The script's "strict" row reproduces the strict mean above (12.21); its medians and 16-84% differ by
  up to 0.1 Gyr because they are taken directly on the 0.2 Gyr grid, not from the histograms.
- **Conclusion:**
  - The isochrones are the dominant difference. With everything else fixed, his isochrones move the best fit from
    14.0 to 12.6 Gyr and the posterior to about 12.2 Gyr. They also fit far better without the map (2.2 vs
    4.3 sigma).
  - Our pipeline with his isochrones gives an age about 1.2 Gyr younger than his published 13.39. His remaining
    choices therefore push his own result older by about that much: his AST error model, 99th-percentile KS,
    Gaussian-process search, weighting and binary fraction. His AST error model alone gave +0.35 Gyr in our
    8-run test.
  - His isochrones pass the calibration stars less often than ours (19% vs 59% with chi2 < 2,
    `../m15_ying_isochrones/`).
