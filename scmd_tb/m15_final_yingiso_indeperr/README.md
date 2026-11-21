# M15 production methodology with Ying's isochrones and independent per-filter errors (2026-11-16)

Purpose: the user asked for Martin's isochrones, the advisor's artificial-star error model (each filter's error drawn
independently from its own (radius, magnitude) distribution), and our methodology for everything else.

| | this run | `../m15_final_indeperr` | production `../m15_final` |
|---|---|---|---|
| isochrones | Ying's (Zenodo 13998637), 8-16 Gyr in 0.2 | our DSEP4, 8-15 in 0.1 | our DSEP4 |
| artificial-star errors | independent per filter | independent per filter | both from one artificial star |
| everything else | production's: zero-point map, E 0-0.12, ridge cut, 1/c-weighted radii, error sign, 95th-percentile KS, grid + KS refinement, bootstrap, strict / shifted / scaled posteriors | same | same |

`setup.json` is production's with `iso`, `iso_format` and `ages` taken from `../m15_yingiso_nozp` and `ast_errors`
"independent". `code/` is identical to production's. `runs.txt` and `refine/setup_override.json` are production's.
`m15yi_statistics.py` is a copy of `../m15_yingiso_nozp`'s (age bins to 16 Gyr).

## Validation (`validate_yingiso_indeperr.py`, `logs/validate_yingiso_indeperr.out`)

1. This code reproduces `../m15_yingiso_nozp`'s best cell (10089, 12.6 Gyr, 15.50, 0.086): KS 0.0134053 vs
   0.0134053 (diff 7e-16).
2. KS at single cells (one sCMD, crn seed):

| cell | production + Ying's isochrones (joint errors) | + independent errors (this setup) |
|---|---|---|
| 10089, 12.6 Gyr, 15.50, 0.086 | 0.01016 | 0.03056 |
| 10357, 14.2 Gyr, 15.58, 0.096 (production's best) | 0.41258 | 0.35728 |
| 10411, 13.2 Gyr, 15.429, 0.080 (replica's best) | 0.04728 | 0.05398 |
| 10336, 13.2 Gyr, 15.443, 0.085 (replica repeat's best) | 0.05031 | 0.05715 |

   With the map and joint errors, Ying's isochrones reach 0.0102 at that cell. That is below production's best (0.0137)
   and inside production's null (0.0077 +/- 0.0026). This one-change run (Ying's isochrones only) has not been done.

## Files

- `run.sh` (driver, `driver.pid`, `run.log`): grid `m15yie_scmd` -> refinement `m15yie_refine` ->
  `../m15_final_yingiso_indeperr_refined/`: bootstrap `boot_m15yie`, statistics, plots, `compare_runs.py`
  (`logs/compare_runs.out`).
- `start_after_indeperr.sh` waits for the `../m15_final_indeperr` driver to finish, because that run holds all 272 CPUs,
  then runs `run.sh`. Queued 2026-11-16.

## Results (finished 2026-11-17 01:57; `logs/compare_runs.out`)

Run times: grid 2026-11-16 19:4x-23:38, refinement to 01:45, bootstrap at the lowest cell (10485, 12.8 Gyr, 15.41,
0.055) to 01:55. No errors.

| run (production methodology) | best fit | sigma | within 3 sigma | shifted age | scaled age |
|---|---|---|---|---|---|
| production (our isochrones, joint errors) | 14.2 Gyr | 2.3 | 256 cells (strict 14.22) | 14.25 | 14.10 |
| + independent errors (`../m15_final_indeperr`) | 14.6 Gyr | 4.2 | 0 | 14.66 | 14.63 |
| + Ying's isochrones + independent errors (this run) | 12.8 Gyr (10485, 15.41, 0.055) | 7.4 | 0 | 13.01 (median 12.56, 11.89-14.49) | 13.31 (median 13.08, 11.84-14.70) |
| reference: Ying's isochrones, joint errors, no map (`../m15_yingiso_nozp`) | 12.6 Gyr | 2.2 | 1,911 (strict 12.21) | 12.31 | 12.39 |

- **Best fit by age** (sigma above this run's null): 8.7 (11.2), 8.0 (12.0), 7.4 (12.8), 8.6 (13.2), 8.2 (14.0),
  7.9 (14.8). The curve is flat, with every age equally poor.
- **Posteriors:** the shifted posterior is bimodal. About 70% of the weight is at 12.4-12.8 Gyr and about 25% at
  14.4-15.0, with a dip at 13.4. Its mean of 13.0 falls between the peaks. The scaled posterior is flat from 11.4 to
  14.6. Neither measures an age, so the closeness of 13.0-13.3 to Ying's 13.39 is not a match.
- **Conclusion:** with Ying's isochrones, independent errors make the fit much worse (2.2 -> 7.4 sigma) and leave the age
  unconstrained. Across all setups, independent errors never fit within 3 sigma.
- **Not yet run:** production with only Ying's isochrones (joint errors, map). It gave KS 0.0102 at 12.6 Gyr in the
  validation, inside production's null.
