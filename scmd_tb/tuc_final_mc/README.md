# tuc_final_mc: 47 Tuc production fit on the Monte Carlo-on DSEP4 isochrones (2026-11-22)

**Purpose.** The production 47 Tuc fit (`../tuc_final` -> `../tuc_final_refined`) used `tuc_iso_d4v2`, where no
DSEP4 run read its Monte Carlo var file (LMonteCarlo unset for [alpha/Fe] != 0). This repeats production with one
change: the isochrones are `tuc_iso_d4v2_mcfix`.
- They come from `tuc_dsep4_batch/run_tuc_dsep4_mcfix.batch`, with LMonteCarlo=.TRUE. for every [alpha/Fe].
- The Monte Carlo line "triple-alpha luminosity reached 1% of total" is present in the logs for every mass of 0.8 Msun
  and up, in every finished run checked.

**Start.** The user asked on 2026-11-22 for this to run after the NGC 6791 W14feh work. `run.sh` waits for
`../../scmd_W14feh_pipeline/run_W14feh.sh` to log "all done" and for its jobs to leave the queue.

**Same as production:**
- `setup.json` is `../tuc_final/setup.json` with only `iso` changed. Everything else is production's: zero-point map,
  PDMF -0.45, binaries 0.009, mu 13.15-13.55, E 0-0.05, keep 12.
- `code/` is a copy of `../tuc_final/code`, identical apart from `__pycache__`.

**Steps** (as `../tuc_yingbox_nozp/run.sh`):
1. Grid. Only runs whose isochrones are built at start time go into `runs.txt`; 994 of 1000 were built on 2026-11-22,
   and the last slow DSEP4 runs may not be included. Age limits come from `age_limits_tucmc.py`.
2. KS refinement, then merge.
3. `../tuc_final_mc_refined/`, with a bootstrap of 10,000 fresh references.
4. Statistics: whole grid, `_D17` (mu 13.15-13.45) and `_D17conv`, figures, and `compare_production.txt` (the age
   posteriors next to production's).

**CPU.** physics 240 minus my other running physics jobs (`NP`), plus 32 on free.

**Log:** `run.log`. **PID:** `driver.pid`.

## Results (driver "all done" 2026-11-24 06:40)

**Grid.** The grid used 996 of 1000 runs; 2-3 DSEP4 runs were still running at the start. The bootstrap was 10,000
fresh references at the best cell: run 10567, 10.2 Gyr, mu 13.50, E 0.040. The null is 0.00537 +/- 0.00174, so the 3σ
line is 0.0106.
- The null is a new, independent file, but matches production's to 4 digits: realisation i uses seed i, and with
  46,539 stars the null depends mostly on the star count.

**Best fit.** KS 0.00912 (2.1σ); production's was 0.00892 (2.0σ).

**Passing.** 78 rows pass 3σ over the whole grid, and 41 inside D17. Strict keeps 22 runs (13 in D17), against 93
(72) in production. The top 5 runs carry 62 % of the strict_D17 weight.

**Age posteriors (Gyr; `compare_production.txt`):**

| | production (Monte Carlo off) | this (Monte Carlo on) |
|---|---|---|
| strict | 10.97 +/- 0.55 (median 10.86, 93 runs) | 11.16 +/- 0.85 (11.01, 22) |
| strict_D17 | 11.40 +/- 0.37 (11.37, 72) | 11.75 +/- 0.63 (11.56, 13) |
| shifted | 11.20 +/- 0.76 (11.11, 320) | 11.15 +/- 0.87 (11.02, 272) |
| shifted_D17 | 11.65 +/- 0.58 (11.55, 248) | 11.72 +/- 0.69 (11.60, 190) |
| scaled | 11.28 +/- 0.88 (11.17, 454) | 11.23 +/- 0.95 (11.10, 520) |
| scaled_D17 | 11.78 +/- 0.70 (11.65, 334) | 11.81 +/- 0.77 (11.68, 382) |

**Other parameters, strict_D17** (production -> this):
- [Fe/H] -0.731 -> -0.749;
- mu 13.416 -> 13.387;
- E 0.028 -> 0.020;
- Y 0.258 -> 0.258;
- mixing length 1.84 -> 1.74;
- [alpha/Fe] 0.40 in both. The grid has only 0.2 and 0.4, and 0.4 fits better.

**Conclusion.**
- Unlike M15, where the age moved 0.3-0.5 Gyr younger, reading the Monte Carlo inputs leaves the 47 Tuc age
  essentially unchanged: within 0.07 Gyr for shifted and scaled, which keep 190-520 runs.
- Strict moves 0.2-0.35 Gyr older (median +0.15-0.19), but rests on 13-22 runs. Its width grows from 0.37 to 0.63 Gyr.
- The varied physics broadens every posterior.

## Comparison with Ying et al. 2025 (ApJ, arXiv:2505.02969, Table 7), 2026-11-24

**Ying et al. 2025, 47 Tuc:**
- 2D KS: 11.76 +/- 0.76 Gyr, [Fe/H] -0.77 +/- 0.03, mu 13.40 +/- 0.05, E 0.02 +/- 0.01.
- Voronoi binning (his adopted method): 11.89 +/- 0.84 Gyr.
- DEBs: 11.36 +/- 0.81 Gyr.

His 47 Tuc priors (Table 4: mu 13.15-13.45, E(V-I) 0-0.05) are the D17 box.

**Ours, strict_D17:** 11.75 +/- 0.63 Gyr, [Fe/H] -0.749, mu 13.387, E 0.020. That is -0.01 Gyr against his 2D KS.
shifted_D17 and scaled_D17 give 11.72 and 11.81.

Production (Monte Carlo off) was 11.40, 0.36 below him. With the Monte Carlo inputs read, the two agree in every
parameter.

**Setup differences that remain:**
- PDMF: ours -0.45, his -0.84.
- [alpha/Fe]: ours fits 0.40 (the grid has 0.2 and 0.4); his fiducial fit gives 0.20-0.24.
- Grid size: ours 1000 runs, his 10,000.
- Error model: ours takes correlated AST errors from one artificial star.
- Weights: his include calibration stars.
