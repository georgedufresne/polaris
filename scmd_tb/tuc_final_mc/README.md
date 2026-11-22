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
