# M15 with Ying's reddening box and no zero-point map (2026-11-10)

Purpose: move the production M15 fit two steps toward M. Ying's setup, as asked on 2026-11-10, and see what
the age does. The production pipeline (`../m15_final/` -> `../m15_final_refined/`) is rerun on the same 500
runs with two changes. Both changes are made together, so this run alone cannot separate their effects (see
"Reading the result").

## The two changes

| | production | this run |
|---|---|---|
| reddening range | grid E 0.050-0.120 (step 0.005), refinement E 0.000-0.120 (step 0.001) | Ying's E(V-I) box 0.08-0.15 applied directly as the F606W-F814W shift: grid 0.080-0.150 (0.005), refinement 0.080-0.150 (0.001) |
| spatial colour zero-point map | on (`zp_correct`, reference F606W 19.3-20.6, F606W share 0.57) | off (`zp_correct` false), as in Ying's paper and KS code |

mu 15.30-15.60 is unchanged: production already used Ying's range. `setup.json` differs from
`../m15_final/setup.json` only in `zp_correct` and `ebv_grid`; `refine/setup_override.json` sets the
refinement E grid.

Not changed (production): ridge cut on data and sCMDs, nearest-AST errors (k 30, crowding metric), radii from
the raw photometry with 1/completeness weights, AST alignment by the zero points zp_v / zp_i, PDMF -1.00,
binaries 0.017, DSEP4 isochrones with MIST BCs, 95th-percentile 2D KS, stage 2 + 3 grid, KS refinement
(8-neighbour descent, not a Bayesian optimiser), 10,000-realisation bootstrap with fresh references at the
lowest-KS cell, strict / shifted / scaled posteriors.

## Files

- `code/`: copy of `../m15_final/code`, unchanged.
- `validate_yingbox.py`, `validate.sbatch` (output `logs/validate_*.out`; it also checks 47 Tuc): this code
  with the production setup reproduces the production best cell; the new setup applies no zero-point map.
- `run.sh` (driver; `driver.pid`, log `run.log`): grid -> collect -> refinement -> merge ->
  `../m15_yingbox_nozp_refined/` -> bootstrap -> statistics -> figures. `../tuc_yingbox_nozp/run.sh` starts
  when this one finishes.
- Outputs: `../m15_yingbox_nozp_refined/posterior_*`, `analysis/` there.

## Validation (logs/validate_11016768.out)

1. Production reproduced: with the production setup, this code gives the m15_final_refined best cell (run
   10357, 14.2 Gyr, mu 15.58, E 0.096) KS 0.013690, as in the file.
2. New setup: zp_correct false, no map applied, refinement E grid 0.080-0.150. 32,089 fit stars (production
   32,085). The same cell without the map has KS 0.0213. That cell was optimised for the corrected
   photometry, so the best cell is expected to move.

## Reading the result

Compare with production's `posterior_strict_DB` (mu 15.3-15.6, E 0.08-0.15 inside production's
E <= 0.12 refinement: 14.24 Gyr). That posterior has the same box with the map on. The difference from it is
roughly the effect of removing the map inside Ying's box. The difference from production's
`posterior_strict` (14.18 Gyr) is the combined effect of both changes.

## Results

Run 2026-11-10 16:58 -> 2026-11-11 04:52 (`run.log`): grid 9.0 h, refinement 2.7 h, bootstrap 9 min.
Outputs in `../m15_yingbox_nozp_refined/`.

- Best fit: run 10245, 14.0 Gyr, mu 15.54, E 0.085, KS 0.01895. The next two are 14.1 and 14.0 Gyr.
  Bootstrap null 0.00762 +/- 0.00261 (3 sigma line 0.01545). The best fit is 4.3 sigma above the null
  (production: 0.0137, 2.3 sigma), so no cell passes 3 sigma and there is no strict posterior.
- Lenient posteriors (production in brackets, then production restricted to Ying's box, `_DB`):

| posterior | this run | production | production `_DB` |
|---|---|---|---|
| strict | none | 14.18 (14.03-14.32) | 14.24 (14.13-14.39) |
| shifted | 14.11 (13.08-14.34) | 14.20 (14.02-14.39) | 14.27 (14.12-14.45) |
| scaled | 13.11 (12.11-14.12) | 14.15 (13.76-14.39) | 14.21 (13.85-14.42) |

- The KS-against-age curve flattens (`analysis/ks_vs_age.png`). Best KS at each age, against production
  restricted to E >= 0.08 (the same box with the map on; the box's upper E edge does not bind):

| age (Gyr) | 11.4 | 12.3 | 13.0 | 13.3 | 14.0 | 14.2 | 14.5 |
|---|---|---|---|---|---|---|---|
| map on (production, E >= 0.08) | 0.0232 | 0.0218 | 0.0232 | 0.0253 | 0.0165 | 0.0137 | 0.0145 |
| map off (this run) | 0.0225 | 0.0202 | 0.0221 | 0.0232 | 0.0190 | 0.0194 | 0.0198 |

  Removing the map makes the 14.0-14.5 Gyr fits worse by 0.003-0.006 and leaves the 11-13 Gyr fits the
  same or slightly better. With the map, 12.3 Gyr is 0.008 worse than the best; without it, only 0.0012.
- The young fits sit on the box edges: mu 15.60 (Ying's upper limit) at 11.4-12.6 Gyr, and E 0.080 (his
  lower limit) at 12.3 Gyr.

Bootstrap from the observed stars themselves (`obs_bootstrap_check.py m15_final/setup.json 50`, output in
`logs/obs_bootstrap_check.out`, 2026-11-11). This uses the production fit stars (32,085 in the KS window),
scored like the fits (95th-percentile 2D KS):

| null | KS mean +/- std | mean + 3 sigma |
|---|---|---|
| ours: model sCMD vs model sCMD at the best cell, artificial stars resampled | 0.0076 +/- 0.0026 | 0.0155 |
| observed stars vs a resample of themselves | 0.0035 +/- 0.0010 | 0.0066 |
| two independent resamples of the observed stars | 0.0059 +/- 0.0019 | 0.0116 |
| two random halves of the observed stars (N/2 each) | 0.0076 +/- 0.0021 | 0.0138 |

The data-only nulls are all narrower than ours. A resample against its own parent shares about 63% of its
stars, so it understates the noise of the data-vs-sCMD comparison, where the two samples are independent.
Two independent resamples are the like-for-like data-only version. Ours is about 30% larger than that,
mostly from the artificial-star resampling, which no data-only bootstrap can include.

Conclusion: in Ying's box, the zero-point map is what separates 14 Gyr from 12-13 Gyr. Without it the best
fit is still at 14.0 Gyr, but the fit is much worse and nearly as good from 11.4 to 14.5 Gyr. The age then
depends on how the posterior is built: 14.1 (shifted) or 13.1 (scaled). This revises the earlier 8-run
estimate (map on vs off, +0.05 Gyr; local grid 12.6-14.8 Gyr, E 0.05-0.10). Which map is used still makes
no difference (`../m15_martinzp/`).
