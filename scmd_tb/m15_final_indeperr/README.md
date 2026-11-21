# M15 production with independent per-filter artificial-star errors (2026-11-16)

Purpose: the advisor recommends treating the F606W and F814W errors as independent, as M. Ying's MakeTestCMD does:
each filter's error is drawn from its own distribution in distance from the centre and magnitude. This is a one-change
test against production (`../m15_final` -> `../m15_final_refined`, 14.2 Gyr).

| | this run | production |
|---|---|---|
| artificial-star errors | `ast_errors` "independent": each filter's completeness test and error drawn independently from its own 10 x 12 (radius, mag) bins, as MakeTestCMD | "nearest": both filters' errors and recovery from one of the 30 nearest artificial stars |
| everything else | production's: DSEP4 isochrones, zero-point map, E 0-0.12, ridge cut, 1/c-weighted radii, error sign, 95th-percentile KS, grid + KS refinement, bootstrap at the lowest cell, strict / shifted / scaled posteriors | |

`setup.json` is `../m15_final/setup.json` with only `ast_errors` changed. `code/` is identical to production's.
`runs.txt` and `refine/setup_override.json` are production's.

Related: `../m15_ast_error_correlation`. There, M15's artificial-star errors in the two filters correlate at about
+0.84, and independent draws predict a main sequence 1.7-1.9x wider than observed.

## Validation (`validate_indeperr.py`, `logs/validate_11083604.out`)

1. This code with production's setup reproduces production's best cell (10357, 14.2 Gyr, mu 15.58, E 0.096): KS
   0.0136905 vs 0.0136905 in its file (diff 4e-15).
2. The same cells with independent errors (one sCMD, crn seed):

| cell | production | independent |
|---|---|---|
| 10357, 14.2 Gyr, 15.58, 0.096 | 0.01369 | 0.03282 |
| 10392, 13.3 Gyr, 15.51, 0.065 | 0.01678 | 0.03245 |

## Files

- `run.sh` (driver, `driver.pid`, `run.log`): grid `m15ie_scmd` -> collect -> refinement `m15ie_refine` -> merge ->
  `../m15_final_indeperr_refined/`: bootstrap `boot_m15ie`, `../m15f_statistics.py`, plots,
  `compare_with_production.py` (`logs/compare_with_production.out`).

Started 2026-11-16 11:52 (grid jobs 11083612 / 11083613).

## Results (finished 2026-11-16 19:44)

Run times: grid 11:52-18:03, refinement to 19:35, bootstrap at the lowest cell (10050, 14.6 Gyr, 15.48, 0.078) to 19:42.
No errors. See `logs/compare_with_production.out` and `../m15_final_yingiso_indeperr/logs/compare_runs.out`.

| | production | this run (independent errors) |
|---|---|---|
| null | 0.00765 +/- 0.00259 | 0.00659 +/- 0.00240 |
| best fit | 0.01369, 2.3 sigma (10357, 14.2 Gyr, 15.58, 0.096) | 0.01675, 4.2 sigma (10050, 14.6 Gyr, 15.48, 0.078) |
| cells within 3 sigma | 256 | 0 (no strict posterior) |
| shifted age | 14.25 (14.02-14.39) | 14.66 (14.47-14.76) |
| scaled age | 14.10 (13.76-14.39) | 14.63 (14.41-14.77) |

**Conclusion:** independent per-filter errors make the fit worse (2.3 -> 4.2 sigma, nothing within 3 sigma) and move the
age older by about 0.4-0.5 Gyr, away from Ying's 13.39. This is the same direction as `../m15_errchain` (+0.35 Gyr).
It agrees with `../m15_ast_error_correlation`: the artificial-star errors of the two filters are correlated, and
independent draws predict a main sequence wider than observed.
