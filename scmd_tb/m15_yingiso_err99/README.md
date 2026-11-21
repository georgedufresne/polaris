# M15: Ying's isochrones + his artificial-star error model + 99th-percentile KS (2026-11-12)

Purpose: take the recreation `../m15_yingiso_nozp` (Ying's published isochrones, no zero-point map, his E box; best
age 12.6 Gyr, posterior about 12.2) one step closer to Ying et al. (2025). This adds his artificial-star error model
and his 99th-percentile 2D KS, and keeps our search (full grid + KS refinement).

| | this run | `../m15_yingiso_nozp` |
|---|---|---|
| isochrones | Ying's (Zenodo 13998637) | same |
| zero-point map; E box | off; mu 15.30-15.60, E 0.08-0.15 applied directly | same |
| artificial-star errors | **his model**: `ast_errors` "independent" (per band, completeness test + error from 10 radial x 12 mag bins, as MakeTestCMD / PhotoPropI), `radial` "observed" (unweighted fit-star radii), `martin_sign` true (observed = true + (input - output), as his Fortran) | ours: errors and completeness copied together from one of the 30 nearest ASTs; 1/c-weighted raw radii; sign observed = true - (input - output) |
| 2D KS | **99th percentile** (his code) | 95th |
| kept ours | grid + KS refinement (not his 50-call Gaussian-process search); ridge cut on data and sCMDs; PDMF -1.00; binaries 0.017 (his 0.013); 32,089-star sCMDs (his 4M); AST placement by photometric zero points (not his SGB matching); bootstrap with fresh references at the lowest cell; strict / shifted / scaled posteriors | same |

The code (`code/`) is an unchanged copy of `../m15_yingiso_nozp/code`. `setup.json` is that run's setup with the four
keys above changed. `runs.txt`, `refine/setup_override.json` and `m15yi_statistics.py` are the same.

## Validation (`validate_err99.py`, `logs/validate_11041604.out`)

1. This code with `../m15_yingiso_nozp_refined/setup.json` replays that run's best cell (10089, 12.6 Gyr, mu 15.50,
   E 0.086): KS 0.0134098 vs 0.0134053 in the file (+4e-6). The file's row comes from the refinement, whose E steps
   are summed (0.08 + k x 0.001), so a rounding-level difference is expected. The code is byte-identical to that
   run's.
2. The same cell with the changes added one at a time (one sCMD, crn seed; the cell was optimised for our model):

| setup | KS |
|---|---|
| `../m15_yingiso_nozp` | 0.0134 |
| + his error model (independent draws, observed radii) | 0.0270 |
| + his error sign | 0.0305 |
| + 99th-percentile KS (this setup) | 0.0458 |

   The 99th-percentile KS is larger by construction, so it is compared with its own bootstrap null.
   `compare_sigma_vs_age.py` plots each run in sigma above its own null.

## Files

- `run.sh` (driver, `driver.pid`, `run.log`): grid `m15ye_scmd` -> collect -> refinement `m15ye_refine` -> merge ->
  `../m15_yingiso_err99_refined/`: bootstrap `boot_m15ye`, statistics, plots, `analysis/sigma_vs_age.png`.

## Results (finished 2026-11-12 22:54)

The run went grid 17:07-21:32 (2,500 units; his error model is about 2x faster per sCMD) -> refinement to 22:46
(0 errors) -> bootstrap to 22:53 -> statistics. Outputs are in `../m15_yingiso_err99_refined/` and
`analysis/sigma_vs_age.png`.

- **Bootstrap null** (run 10027, 12.8 Gyr, mu 15.47, E 0.081; 99th-percentile KS): 0.00759 +/- 0.00233, 3 sigma
  line 0.01459. Its 95th-percentile column is 0.0062, against 0.0077 for the null with our error model, so his
  error model gives a narrower null.
- **Fit quality:**
  - Best KS 0.02350 (run 10027, 12.8 Gyr, mu 15.47, E 0.081), 6.8 sigma above the null. No cell is within
    3 sigma, so there is no strict posterior.
  - Next: 10189 at 13.2 Gyr (0.02432) and 10213 at 12.6 Gyr (0.02433).
  - At the 95th percentile the best fits are 0.0213, against 0.0134 with our error model. Most of the loss of fit
    therefore comes from his error model, as in `../m15_errchain` (his model fitted worst there too).
- **Best KS by age (sigma above own null):**

| age (Gyr) | 11.2 | 11.6 | 12.0 | 12.4 | 12.8 | 13.2 | 13.6 | 14.0 | 14.4 | 14.8 |
|---|---|---|---|---|---|---|---|---|---|---|
| his isochrones, our errors, p95 (`../m15_yingiso_nozp`) | 2.5 | 2.3 | 2.2 | 2.4 | 2.4 | 2.3 | 3.0 | 3.8 | 3.9 | 4.6 |
| his isochrones, his errors, p99 (this run) | 8.0 | 8.2 | 7.3 | 7.2 | 6.8 | 7.2 | 9.0 | 10.7 | 10.5 | 9.5 |

- **Posteriors (age in Gyr; `logs/calib_reweight.out` for the calibration weights):**

| | shifted | scaled |
|---|---|---|
| KS weights | 12.42 (median 12.35, 16-84% 11.84-12.82) | 12.57 (median 12.37, 11.61-13.26) |
| times calibration weights (his isochrones) | 12.33 (11.80-12.80) | 12.22 (11.40-13.00) |

- **No-cut weighting** (`nocut_posterior.py`, `logs/nocut_posterior.out`, 2026-11-14). This is the rule of Ying et
  al. (2024, 2025) for the 2D KS: each isochrone gets "the probability" of its KS value in the null, with no threshold.
  The papers give no formula, so three readings were tried (KDE density, fitted normal, empirical survival fraction),
  per (run, age) cell and per run:
  - `../m15_yingiso_nozp` (our errors): 12.23-12.27 Gyr, or 12.12-12.17 with calibration weights. This equals its
    strict posterior (12.21).
  - This run: every fit is outside the null (null maximum 0.0203, best fit 0.0235).
    - The KDE reading puts all the weight on the best cell: 12.8 Gyr.
    - The normal reading gives 12.7-12.8 (12.5 with calibration weights), carried by 2-4 runs.
    - The survival reading gives zero weight everywhere.
  - So without a cut no age goes above the best fits' 12.8. Only the KS values themselves, not the weighting, can
    move the age to 13.4.
- **Conclusion:**
  - His error model plus the 99th-percentile KS move the age older by about 0.2 Gyr: best fit 12.6 -> 12.8 Gyr,
    posteriors about 12.2-12.4 -> 12.4-12.6. The fit gets much worse (2.2 -> 6.8 sigma).
  - The best-fit distance and reddening now match Ying's published values (DM 15.47-15.50 vs his 15.49 +/- 0.03;
    E 0.081-0.083 vs his 0.08 +/- 0.01).
  - The age is still about 0.8-1.0 Gyr below his 13.39.
  - What differs still: his Gaussian-process search, his weighting (null density at every isochrone, no cut), 4M-star
    sCMDs and a 1M-star null, his 0.08-mag cut around each isochrone (ours: the ridge cut), binaries 0.013, and
    SGB-matched AST placement. Voronoi is his main M15 method (13.23).
