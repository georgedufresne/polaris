# M15 production with PHOENIX bolometric corrections: the one-change colour-table test (2026-11-17)

Purpose: `../m15_isochrone_bc_test` found that Ying's isochrones differ from ours in two ways. First, the colour tables:
his PHOENIX corrections are 0.035-0.056 mag redder than our MIST ones and explain his ~0.03 redder isochrones. Second,
the stellar models: the per-run turnoff and SGB structure differs even with the same tables. This run converts our own
models with PHOENIX and leaves the rest of production unchanged. Together with `../m15_final_yingiso` (his isochrones),
it splits the isochrone effect on the age into a colour-table part and a model part.

| | this run | production `../m15_final` | `../m15_final_yingiso` |
|---|---|---|---|
| stellar models | our DSEP4 | our DSEP4 | Ying's |
| colour tables | Dartmouth PHOENIX ACS/WFC ([Fe/H] -2.25, [alpha/Fe] +0.4) | MIST ACS/WFC ([Fe/H] -2.25) | his (PHOENIX, nearest [alpha/Fe]) |
| everything else | production's | | production's |

- `make_phoenix_bc.py` (`logs/make_phoenix_bc.out`) writes `inputs/bc_phoenix_acs_wfc_fehm225_ap4.npz`: the mean of the
  [Fe/H] -2.0 and -2.5 tables on their 789 shared (Teff, log g) points, with bc = -BCv + (V - X) (in these tables
  BCv = M_V - M_bol). On runs 10357, 10392 and 10089 it agrees with `../m15_isochrone_bc_test`'s conversion to 1.7 mmag
  in F606W and 1.1 mmag in colour. The effect under test is 35-56 mmag.
- `setup.json` is production's with `bc` changed. `code/`, `runs.txt` and `refine/setup_override.json` are production's.
- Checks (`logs/validate_cells.out`, one sCMD per cell):

| cell | production (MIST) | PHOENIX |
|---|---|---|
| 10357, 14.2 Gyr, 15.58, E 0.096 (production's best) | 0.01369 | 0.31431 |
| same at E 0.060 | 0.21648 | 0.02203 |
| 10392, 13.3 Gyr, 15.51, E 0.030 | 0.20404 | 0.02482 |

  PHOENIX shifts the fit to about 0.036 lower E, inside the 0-0.12 range.
- `run.sh` (jobs `m15px_*`, `../m15_final_phxbc_refined/`, `compare_runs.py` -> `logs/compare_runs.out`).
- `start_after_yingiso.sh` waits for the `../m15_final_yingiso` driver, then runs `run.sh`.
