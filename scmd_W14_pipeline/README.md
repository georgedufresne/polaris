# NGC 6791 sCMD leg with W v1.4 isochrones (2026-11-18 ->)

Purpose: the NGC 6791 sCMD fit with the updated W emulator (CONF1DENCE-W fcfac3c, W v1.4). The fit is compared with the
other legs: OG, DSEP3 with the physics on ("New", `../analysis_noco`), the old emulator, and W v1.2 (October).
Isochrones: `../emul_W14` (README there).

## Pipeline

- The scripts are copies of `../scmd_emul_pipeline` (the W leg). Only the folder names `Isochrones_W -> Isochrones_W14`,
  `sCMD_W -> sCMD_W14`, `2d_ks_W -> 2d_ks_W14`, `plots_for_2d_W -> plots_for_2d_W14` and `analysis_W -> analysis_W14`
  change; `diff` shows nothing else.
- `run_scmd_W14.batch` is `run_scmd_W_physics.batch` with the W14 folder and a run list.
- `merge_W14.py` and `prefilter_W14.py` keep the fixed KS cut 0.0468, the 3 sigma line of the original bootstrap null,
  used for every leg.
- `run_W14.sh` (driver): isochrones -> sCMD array (40 at once on physics) -> merge -> prefilter ->
  `../scmd_noco_pipeline/analyze_scmd_only.py` -> `scmd_5way_compare.py` (-> `../scmd_5way_comparison.txt`).
- **KTTAU = 1 excluded.** W v1.4 is not trained on KTTAU = 1. Its KTTAU = 1 sets fit only at 10-11.5 Gyr while DSEP3
  fits the same sets at 8.5-9 Gyr.
  - 278 queued KTTAU = 1 tasks were cancelled on 2026-11-19 (`kt1_cancelled.txt`). `finish_watch.sh` releases the driver
    when the other 722 tasks are done.
  - `scmd_5way_compare.py` and `age_breakdown.py` restrict EVERY leg to KTTAU 0/5. The KTTAU of each run is the same in
    every leg's `initial_vars.csv` and in the var files.

## Interim results (2026-11-21 14:00; 403 W v1.4 runs with KS results)

The same 290 runs (KTTAU 0: 145, 5: 145) are fitted in all three legs; each leg applies its own KS cut:

| | OG | DSEP3 physics on | W v1.4 |
|---|---|---|---|
| runs passing the cut | 168 (58%) | 171 (59%) | 166 (57%) |
| age, mean +/- sd (Gyr) | 8.51 +/- 0.62 | 8.51 +/- 0.65 | 8.36 +/- 0.65 |
| median (16-84%) | 8.49 (7.85-9.08) | 8.47 (7.79-9.16) | 8.31 (7.65-8.99) |
| the 133 runs kept by all three | 8.48 +/- 0.61 | 8.50 +/- 0.64 | 8.34 +/- 0.63 |

- **Age differences** (runs resampled together, 1000 resamples):
  - W v1.4 - DSEP3: -0.143 +/- 0.030 Gyr (z -4.7).
  - W v1.4 - OG: -0.143 +/- 0.043 (z -3.3).
  - DSEP3 - OG: +0.001 +/- 0.041 (z 0.0).
  - On 2026-11-20 (175 runs) these were -0.13 +/- 0.04, -0.12 +/- 0.05 and +0.02 +/- 0.05.
- **Hung tasks:** 24 tasks had logged `done` but stayed RUNNING (the Slurm spank-auks bug); they were cancelled on
  2026-11-21. `reap_zombies.sh` (pid in `reap_zombies.pid`, log `zombies_cancelled.txt`) now cancels such tasks every
  10 min, and only after their `done` line is in the log.
- **Why W v1.4 is younger** (same runs, 8.5 Gyr, 100 runs, median (16-84%)):
  - Turnoff log L agrees (+0.000 dex; -0.016..+0.013).
  - Turnoff log Teff is -0.0021 dex (-0.0032..-0.0007), about 27 K cooler. The upper MS (log L 0.2) is -0.0020 dex.
  - A redder isochrone at fixed age looks older, so the fit takes a slightly younger age (and lower E(B-V):
    0.177 vs 0.181).
  - It is not DSEP3 vs DSEP4: the 47 Tuc comparison with the physics on agrees to 2e-4 dex at the turnoff,
    `../scmd_tb/tuc_iso_compare/README_mcfix.md`.
  - It is not the conductive-opacity convention: at most 6e-4 dex in the fit window, `../emul_W14/README.md`.
  - It is an emulator temperature bias of ~0.002 dex at the turnoff.

`age_breakdown.py` (`age_breakdown.txt`, `.png`) gives the posterior shapes, per-set differences and fit quality,
KTTAU 0/5 only.
