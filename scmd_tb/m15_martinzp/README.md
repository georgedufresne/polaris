# M15 with Martin's zero point (2026-11-09)

Purpose: check that our colour zero-point logic does not change the M15 age compared with M. Ying's. The
production pipeline (`../m15_final/` -> `../m15_final_refined/`) is rerun on the same 500 runs with one
change: Martin's zero-point correction in place of ours. Everything else, the statistics included, is
production.

## The one change

| | ours (production) | Martin (this test) |
|---|---|---|
| stars used | ridge-selected upper MS, F606W 19.3-20.6 | every star with F606W within 3.5 mag of the SGB (18.64) |
| offset of a star | colour - ridge colour, minus the AST-predicted crowding shift | colour - fiducial colour (cubic spline of `fiducial_lines.csv`) |
| correction | median of the 50 nearest reference stars, self excluded | median of the 100 nearest stars, self included |
| applied to | colour; F606W takes 0.57 of it | colour only, F606W fixed |

Source: `Files_from_Martin/github_Age-of-GCs/M15_data/M15_ZPCR.ipynb`. Martin's own corrected files
(`M15_fitstars_ZPCR.dat`, `M15_fitstars_2dKS_ZPCR.dat`) are not used: every one of their stars matches
`M15_2023.dat`, none matches our `M15.dat`. The notebook's recipe is applied to our `M15.dat` instead.
Sky positions (median RA/Dec of each catalogue against the Harris 2010 centres, 2026-11-10): our `M15.dat`
is 0.2 arcmin from M15 (NGC 7078) and 13 deg from M2; `M15_2023.dat` is 0.1 arcmin from M2 (NGC 7089) and
13 deg from M15. Martin's ZPCR files are M2 stars labelled M15. Matched by position and F606W against both
catalogues (2026-11-10): `M15_fitstars_ZPCR.dat`, `M15_fitstars_2dKS_ZPCR.dat` and `M15_fitstars_with_bins.dat`
are M2; `M15_fitstars.dat`, `M15_fitstars_cmd.dat` and `M15_fitstars_2dKS.dat` are M15. His KS code reads
`M15_fitstars.dat`, and the paper's methods mention no zero-point correction.

Code: `code/` is a copy of `../m15_final/code` plus `ast.martin_zero_point` and `Setup.zp_method`
("ours" | "martin", default "ours", which reproduces production). `setup.json` differs from
`../m15_final/setup.json` only in `zp_method` "martin", `zp_k` 100, `zp_mag_frac` 0, `zp_martin_dmag` 3.5.

Not changed (production): ridge cut on data and sCMDs, nearest-AST errors (k 30, crowding metric), radii
from the raw photometry with 1/completeness weights, PDMF -1.00, binaries 0.017, DSEP4 isochrones with MIST
BCs, stage 2 + 3 grid, KS refinement (E step 0.001, E 0.000-0.120), 10,000-realisation bootstrap at the
lowest-KS cell, strict / shifted / scaled posteriors, the prior boxes.

## Files

- `validate_martinzp.py` (output in `logs/validate_*.out`): (1) the function reproduces Martin's own
  corrected colours from `check_m2/M15_2023.dat` (copied from the earlier session's scratchpad) and the
  notebook's neighbour loop; (2) with "ours" the code copy reproduces the production best cell; (3) our map
  against Martin's on our data. `analysis/zp_maps.npz`: both maps for every star of `M15.dat`.
- `run.sh` (driver, `driver.pid`, log `run.log`): grid -> collect -> refinement -> merge ->
  `../m15_martinzp_refined/` -> bootstrap -> statistics -> figures.
- Outputs: `../m15_martinzp_refined/posterior_*`, `analysis/` there.
- `ying2025_match_prompt.md` (2026-11-10): a self-contained brief for another AI model on how to match
  Ying et al. (2025) for M15. It summarises this test, `../m15_errchain/`, `../tuc_errchain/`,
  `../m15_final_refined/` and the paper's sections 1-3.

## Validation (logs/validate_11000243.out)

1. Martin's own output reproduced: on `check_m2/M15_2023.dat` with V_SGB 18.85 and
   `check_m2/fiducial_2023_M2.csv`, `martin_zero_point` gives his `M15_fitstars_ZPCR.dat` colours for all
   91,945 stars to 2e-16; the notebook's neighbour loop agrees exactly (300 stars of our `M15.dat`).
2. Production reproduced: with zp_method "ours" this code gives the m15_final_refined best cell (run 10357,
   14.2 Gyr, mu 15.58, E 0.096) KS 0.013690, as in the file.
3. Our data, stars in the window (F606W 16.64-20.64): ours rms 0.0088 mag, Martin rms 0.0086; correlation
   0.70; Martin - ours rms 0.0068, median +0.0008. Fit stars 32,085 (ours) vs 32,116 (Martin). The
   production best cell with Martin's zero point: KS 0.0170 (0.0137 with ours); the best cell may move.

## Results

Run 2026-11-09 17:09 -> 2026-11-10 01:48 (`run.log`): grid 6.8 h, refinement 1.7 h, bootstrap 7 min.
Outputs in `../m15_martinzp_refined/`. Strict posteriors (KS <= null mean + 3 sigma), medians and 16-84%:

| | production (our zero point) | Martin's zero point |
|---|---|---|
| best KS (null 3 sigma line) | 0.0137 (0.0154), 2.3 sigma | 0.0129 (0.0153), 2.1 sigma, p 0.044 |
| cells / runs inside 3 sigma | 149 / 136 | 291 / 248 |
| age, all E | 14.18 (14.03-14.32) | 14.21 (14.06-14.38) |
| age, box mu 15.3-15.6, E 0.062-0.117 | 14.20 (14.05-14.33) | 14.24 (14.08-14.40) |
| [Fe/H], all E | -2.31 (-2.39, -2.21) | -2.31 (-2.40, -2.22) |
| E(606-814), all E | 0.076 (0.065-0.087) | 0.071 (0.058-0.081) |
| mu, all E | 15.515 | 15.506 |

Paired over the 500 runs: each run's best age is unchanged (median difference 0.00 Gyr, within 0.1 Gyr for
99% of runs); Martin's zero point lowers each run's best KS slightly (median -0.0006, lower in 85% of runs),
so more runs pass the same 3 sigma line; E lower by 0.001 (median), mu unchanged. The bootstrap nulls are the
same (0.00762 +/- 0.00256 vs 0.00765 +/- 0.00259).

Conclusion: our zero-point logic does not change the M15 age (+0.03-0.04 Gyr with Martin's, a quarter of
the posterior width). The gap to Ying et al. (13.2-13.4 Gyr) does not come from which zero-point recipe is
used. Having no map at all does matter. In `../m15_yingbox_nozp/` (2026-11-11) the fit gets much worse, and
12-13 Gyr becomes nearly as good as 14 Gyr.
