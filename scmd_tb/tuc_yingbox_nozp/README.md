# 47 Tuc with no zero-point map (2026-11-10)

Purpose: the 47 Tuc counterpart of `../m15_yingbox_nozp/`. The production pipeline (`../tuc_final/` ->
`../tuc_final_refined/`) is rerun on the same 1000 runs without the spatial colour zero-point map, as in
M. Ying's paper and code.

## The change

| | production | this run |
|---|---|---|
| spatial colour zero-point map | on (`zp_correct`, reference F606W 17.9-19.0, F606W share 0.08) | off (`zp_correct` false) |

The M15 run's other change, Ying's reddening box applied directly, has no 47 Tuc counterpart. Production's
E range, 0.000-0.050 applied directly as the F606W-F814W shift, already is Ying's 47 Tuc box. His mu range
13.15-13.45 is applied in the `_D17` posteriors, as it was for production; production's "with priors" age,
11.37 Gyr, is `posterior_strict_D17`. `setup.json` differs from `../tuc_final/setup.json` only in
`zp_correct`. `refine/setup_override.json` is production's (E 0.000-0.050 in 0.001).

Not changed: everything else in production, the statistics included (see `../m15_yingbox_nozp/README.md`).
PDMF -0.45, binaries 0.009.

## Files

- `code/`: copy of `../tuc_final/code`, unchanged.
- `validate_yingbox.py`: run by `../m15_yingbox_nozp/validate.sbatch` (output in that folder's `logs/`).
- `run.sh` (driver; `driver.pid`, log `run.log`): waits for `../m15_yingbox_nozp/run.sh`, then grid ->
  collect -> refinement -> merge -> `../tuc_yingbox_nozp_refined/` -> bootstrap -> statistics (whole grid,
  `_D17`, `_D17conv`) -> figures.

## Validation (../m15_yingbox_nozp/logs/validate_11016768.out)

1. Production reproduced: with the production setup, this code gives the tuc_final_refined best cell (run
   10239, 10.4 Gyr, mu 13.49, E 0.042) KS 0.008918, as in the file.
2. New setup: no map applied. 47,413 fit stars (production 47,418). The same cell without the map has KS
   0.0173. That cell was optimised for the corrected photometry, so the best cell is expected to move.

## Results

**Cancelled by the user on 2026-11-11 at about 14:00**, during the grid stage. Jobs 11020428 and 11020429 were
cancelled and the driver stopped. 4,027 units had finished: 35,077 (run, age) pairs from 585 of the 1,000 runs.
The partial units are kept in `units/`. There was no refinement and no bootstrap of its own.

Partial check (`logs/partial_check.out`), against production's 47 Tuc null (0.0054 +/- 0.0017, 3 sigma line
0.0106):

| | best KS | sigma above null | (run, age) pairs within 3 sigma |
|---|---|---|---|
| no map, whole grid (E 0-0.05, mu 13.15-13.55) | 0.0154 (11.8 Gyr) | 5.8 | 0 |
| no map, mu 13.15-13.45 | 0.0155 (12.6 Gyr) | 5.8 | 0 |
| production, grid stage, same (run, age) pairs | 0.0091 (10.6 Gyr) | 2.2 | 75 |

- **No fit reaches the null without the map.**
- **The KS-age curve is flat without it.** The best KS is 0.0155-0.0169 at every age from 10.2 to 12.9 Gyr,
  against 0.0092-0.0116 in production on the same pairs.

This is the same pattern as M15 (`../m15_yingbox_nozp`). Refinement lowered production's best by about 0.0002,
so it would not change this.
