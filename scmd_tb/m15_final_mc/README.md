# M15 production on the DSEP4 grid rebuilt with the Monte Carlo read fixed (2026-11-18)

Purpose: `../m15_dsep4_montecarlo_check` found that in the original DSEP4 M15 grid only the [alpha/Fe] = 0 control
namelist set `LMonteCarlo`. 484 of 500 runs therefore never read their var file: no varied surface boundary condition,
diffusion, overshoot, nuclear rates or opacities. This reruns the 500 runs with the read on, rebuilds the isochrones,
and repeats production (grid + KS refinement, bootstrap, strict / shifted / scaled posteriors) on them.

## Model grid (M15/dsep4)

- `env/nml_mc/`: a copy of `env/nml` with ` LMonteCarlo=.TRUE.` added before `$END` in the five alpha-enhanced
  `*_3_DATA.nml` files (`diff -r` shows only those lines). The original `env/nml` is untouched.
- `run_m15_dsep4mc.batch`: a copy of `run_m15_dsep4.batch`. `NML` points at `env/nml_mc`, tracks go to
  `M15/tracks_dsep4mc`, isochrones to `M15/iso_dsep4mc`, logs to `logs_mcfix/`. The var files are unchanged
  (ALPHAC = 0 as before; Ying's are 0.03-0.19).
- One-track test (`test_mcfix/`, run 10357 [alpha/Fe] +0.2, 0.78 Msun): the track ends with "triple-alpha luminosity
  reached 1% of total". DSEP prints that only inside `if (lmontecarlo)`, so the var file is now read. Exit 0, 6.6 min.
- Submitted 2026-11-18 13:30 (`jobids_mcfix.txt`): physics array 11103596 (10000-10467, at most 200 at once) and free
  array 11103597 (10468-10499). The original runs that read their var file took about 5.5 h each, so about 14 h in total.

## Fit (this folder)

- `setup.json` is production's with only `iso` changed, to `M15/iso_dsep4mc/{run}/isos/isochrones_example.txt`.
  `code/` and `refine/setup_override.json` are production's.
- `start_after_dsep.sh` (`run.log`):
  - waits until all 500 DSEP tasks log their end, and cancels tasks that hang at exit;
  - writes `runs.txt` (the runs with isochrones) and `m15_age_limits_d4mc.json` (`age_limits_m15mc.py`, a copy of the
    original reading the new folder);
  - then runs `run.sh`: grid `m15mc_scmd` -> refinement -> `../m15_final_mc_refined/` (bootstrap at the lowest cell,
    `../m15f_statistics.py`, plots) -> `compare_runs.py` (`logs/compare_runs.out`, against production and the other
    runs).

**CPU share (2026-11-18, before `run.sh` started):** `run.sh` uses 200 physics + 32 free workers for the grid,
refinement and bootstrap (bootstrap `--ntasks 232`, free `--task-offset 200`, waits for 232 parts) instead of 240 + 32.
That leaves 40 physics CPUs for the NGC 6791 W v1.4 run within the 240 + 32 limit. The original is kept as
`run.sh.272cpu_orig`. The bootstrap splits realisations as `range(task, n, ntasks)`, so 232 tasks cover all 10,000.

## Results (finished 2026-11-19 17:03)

- **Model grid:** the DSEP4 rerun finished 2026-11-19 05:25 (500 of 500 runs, 500 isochrone sets; age limits: 195 runs
  usable to 16 Gyr, the lowest limit 15.0).
- **Fit:** grid 05:27-14:59, refinement to 16:55 (0 errors), bootstrap at the lowest cell (10043, 13.2 Gyr, mu 15.44,
  E 0.054) to 17:02, then statistics. See `logs/compare_runs.out`.
- **Null:** 0.00761 +/- 0.00258 (3 sigma line 0.01536).
- **Best fit:** KS 0.01340 at run 10043, 13.2 Gyr, mu 15.440, E 0.054. That is 2.2 sigma, inside the null. 122 cells are
  within 3 sigma.
- **Best fit at each age** (sigma above own null):

| age (Gyr) | 11.6 | 12.0 | 12.4 | 12.8 | 13.2 | 13.6 | 14.0 | 14.4 | 14.8 |
|---|---|---|---|---|---|---|---|---|---|
| production (Monte Carlo inputs ignored in 484 runs) | 6.5 | 7.1 | 4.8 | 4.1 | 3.3 | 2.8 | 2.5 | 2.7 | 4.2 |
| this run (Monte Carlo inputs read) | 4.3 | 3.9 | 3.5 | 2.6 | 2.2 | 2.6 | 2.6 | 2.5 | 3.2 |

- **Posteriors** (mean +/- sd; median; 16-84%):

| | strict | shifted | scaled | production strict |
|---|---|---|---|---|
| age (Gyr) | 13.75 +/- 0.43 (13.78; 13.20-14.10) | 13.93 +/- 0.61 (14.00; 13.33-14.41) | 13.78 +/- 0.75 (13.90; 12.92-14.43) | 14.22 +/- 0.18 (14.18; 14.03-14.32) |
| [Fe/H] | -2.27 +/- 0.08 (-2.31) | -2.25 +/- 0.10 | -2.25 +/- 0.10 | -2.29 +/- 0.10 |
| DM (F606W) | 15.477 +/- 0.058 | 15.471 +/- 0.054 | 15.480 +/- 0.062 | 15.520 +/- 0.035 |
| E(V-I) | 0.062 +/- 0.020 | 0.057 +/- 0.020 | 0.060 +/- 0.021 | 0.077 +/- 0.012 |
| Y | 0.247 +/- 0.002 | 0.247 +/- 0.002 | 0.247 +/- 0.002 | 0.248 +/- 0.001 |
| [alpha/Fe] | 0.32 +/- 0.12 (0.40) | 0.31 +/- 0.11 | 0.31 +/- 0.11 | 0.28 +/- 0.11 (0.20) |
| mixing length | 1.56 +/- 0.44 (1.34) | 1.57 +/- 0.39 | 1.59 +/- 0.38 | 1.98 +/- 0.34 |
| runs | 43 | 387 | 446 | 136 |

- **Conclusion:**
  - Reading the Monte Carlo inputs moves the M15 age younger by 0.3-0.5 Gyr: strict 14.22 -> 13.75, shifted
    14.25 -> 13.93, scaled 14.10 -> 13.78.
  - The uncertainty roughly doubles, because the physics now varies between runs.
  - The fit improves at 12-13.6 Gyr, and the best fit moves to 13.2 Gyr at 2.2 sigma.
  - Ying et al.'s 13.39 +/- 0.39 is now about 0.6 sigma from the strict posterior (combined errors).

## Comparison with Ying et al. 2025 (ApJ, arXiv:2505.02969, Table 7), 2026-11-24

**Ying et al. 2025, M15:**
- 2D KS: 13.39 +/- 0.39 Gyr, [Fe/H] -2.30, mu 15.49, E 0.08.
- Voronoi binning (his adopted M15 method): 13.23 +/- 0.51 Gyr.

**Ours (strict):** 13.75 +/- 0.43 Gyr, so +0.36 Gyr, 0.6σ combined against his 2D KS.

**Inside his prior box** (mu 15.3-15.6, E(V-I) 0.08-0.15; `EVI_MIN=0.08 EVI_MAX=0.15 MU_MIN=15.3 MU_MAX=15.6
POST_SUFFIX=_yingbox ../m15f_statistics.py`, written to `../m15_final_mc_refined/posterior_*_yingbox/`). The grid's E
reaches 0.12.
- Best fit 0.01414 (2.5σ).
- strict 14.21 +/- 0.18 Gyr (12 runs), shifted 13.89 +/- 0.84, scaled 13.32 +/- 0.94.
- E 0.088-0.093, mu 15.53-15.55.

His reddening floor does not explain the difference. As before (`../m15_ying_isochrones`, `../m15_yingiso_nozp`), the
isochrones are the main lever.
