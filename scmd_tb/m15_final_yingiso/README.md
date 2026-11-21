# M15 production with Ying's isochrones: the one-change isochrone test (2026-11-17)

Purpose: compare our isochrones with Ying's, everything else production. Production (`../m15_final`) gives 14.2 Gyr.
Ying's isochrones gave 12.2-12.6 Gyr in our pipeline without the map and in his E box (`../m15_yingiso_nozp`). His own
code gives 13.2 (`../m15_martin_replica`). This run uses production's methodology exactly, map included, and changes
only the isochrones.

| | this run | production |
|---|---|---|
| isochrones | Ying's (Zenodo 13998637), 8-16 Gyr in 0.2 | our DSEP4 + MIST ACS/WFC BCs, 8-15 in 0.1 |
| everything else | production's: zero-point map, E 0-0.12, ridge cut, errors of one artificial star for both filters, 1/c-weighted radii, 95th-percentile KS, grid + KS refinement, bootstrap, strict / shifted / scaled posteriors | |

`setup.json` is production's with `iso`, `iso_format` and `ages` from `../m15_yingiso_nozp`. It equals
`../m15_final_yingiso_indeperr/setup.json` except `ast_errors`. `code/` is identical to production's. `runs.txt` and
`refine/setup_override.json` are production's. `m15yi_statistics.py` has age bins to 16 Gyr.

Check (`logs/validate_cell.out`): cell 10089, 12.6 Gyr, 15.50, 0.086 gives KS 0.01016, the same as in
`../m15_final_yingiso_indeperr`'s validation with joint errors. That is below production's best (0.01369) and inside
production's null (0.0077 +/- 0.0026). The code reproduces production's and `../m15_yingiso_nozp`'s best cells exactly
(see those validations).

Files: `run.sh` (driver, `driver.pid`, `run.log`; jobs `m15yj_scmd`, `m15yj_refine`, `boot_m15yj`), `compare_runs.py`
(`logs/compare_runs.out`).
