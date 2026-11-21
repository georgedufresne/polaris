# M15 error-model chain (2026-11-04)

Purpose: test how the AST error model and a spatial photometric zero-point map change the M15 sCMD fit,
one change at a time, starting from Martin's model. Same normalised-perpendicular ridge cut for the data
and every sCMD (unchanged from production).

Inputs: `setup.json` (from `m15_grid_d4`), PDMF slope -1.00 and binary fraction 0.017
(`../literature/scmd_cluster_inputs.csv`), 8 runs (`runs.txt`), local grid ages 12.6-14.8,
mu 15.44-15.58 (step 0.035), E(F606W-F814W) 0.050-0.100. Code: `code/` (default settings reproduce
production exactly).

Variants (`chain.py`): V0_martin (per-band errors, unweighted radii), V1_joint (+ joint errors and
completeness), V2_weighted (+ 1/c radii, production), V3_rawradial (+ radii from before the ridge cut),
V4_floor (c floor 0.02), V5_zpcorr (+ zero-point map on colour, reference 19.3-20.6), V6_zpsplit
(+ zero-point map on both bands, F606W share 0.57 from `diag/diag3.py`), and sensitivity runs S_k15,
S_k60, S_mag_half, S_mag_dbl, S_split. `fits_superseded/` holds the scale runs that changed both scales
together (k-NN depends only on their ratio, so they duplicated V3).

Outputs: `fits/<variant>/`, `nulls/<variant>.csv` (240 bootstrap realisations), `chain_summary.csv`,
`chain_summary.png` (`summary.py`), `diag/` (ridge-cut efficiency, completeness floor, neighbourhood
reach, zero points, closure test, reddening-vs-zero-point test, width by radius).

## Results

| variant | best KS | sigma above its null | best cell (run, age, mu, E) |
|---|---|---|---|
| V0 Martin | 0.0188 | 6.4 | 10195, 14.4, 15.510, 0.070 |
| V1 joint | 0.0180 | 4.6 | 10392, 13.3, 15.510, 0.065 |
| V2 production | 0.0169 | 3.7 | 10392, 13.3, 15.510, 0.065 |
| V3 raw radii | 0.0169 | 3.5 | 10392, 13.3, 15.510, 0.065 |
| V5 zero-point map, colour | 0.0157 | 3.1 | 10142, 14.1, 15.545, 0.080 |
| V6 zero-point map, both bands | 0.0153 | 3.0 | 10392, 13.6, 15.475, 0.060 |

- Martin's model fits worst and shifts the per-run best age by +0.35 Gyr (median).
- k, magnitude scale and short/long split change the best KS by at most 0.0005; no age change.
- The spatial pattern is a zero point (F606W share about 0.57), not reddening (share 2.8 excluded).
- Outer colour width (`diag/diag4.py`): missing spread beyond 700 px drops from about 0.012 to about
  0.007 mag with the map; the remaining ~0.007 mag is still unexplained.

Settings carried to `../m15_final/setup.json`: PDMF -1.00, binaries 0.017, radial_source raw,
zp_correct with zp_ref 19.3-20.6 and zp_mag_frac 0.57.
