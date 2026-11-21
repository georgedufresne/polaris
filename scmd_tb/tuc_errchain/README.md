# 47 Tuc error-model chain (2026-11-04)

Purpose: test how the AST error model and a spatial photometric zero-point map change the 47 Tuc sCMD
fit, one change at a time, starting from Martin's model. Same normalised-perpendicular ridge cut for the
data and every sCMD (unchanged from production). Companion to `../m15_errchain/`.

Inputs: `setup.json` (from `tuc_grid_d4`), PDMF slope -0.45 and binary fraction 0.009
(`../literature/scmd_cluster_inputs.csv`), the 8 lowest-KS runs of the production grid (`top8.txt` in
`chain.py`), local grid ages 11.0-13.6, mu 13.300-13.475 (step 0.035), E(F606W-F814W) 0.000-0.050.
Code: `code/` (copy of the M15 experiment code; default settings reproduce production exactly).

Variants (`chain.py`):
- V0_martin: per-band errors and completeness product, unweighted radii
- V2_weighted: joint (dV, dI) errors, 1/c-weighted radii (production)
- V3_rawradial: radii from the photometry before the ridge cut
- V6_zpsplit: + zero-point map, reference upper MS F606W 17.9-19.0, F606W share 0.08 (`TUC_ZPK`; from
  `diag/diag_generic.py` with `diag/tuc.json`, per-star SGB-vs-local-MS slope 2.54 +/- 0.2 against 2.44
  for an F814W-only zero point and 6.94 for reddening)

Outputs: `fits/<variant>/`, `nulls/<variant>.csv` (240 bootstrap realisations at each variant's best
cell), `chain_summary.csv`, `chain_summary.png` (`summary.py`), `diag/`.

## Results

| variant | best KS | null mean + 3 sigma | sigma above null | best cell (run, age, mu, E) |
|---|---|---|---|---|
| V0 Martin | 0.0172 | 0.0105 | 6.7 | 10745, 12.0, 13.405, 0.020 |
| V2 production | 0.0153 | 0.0111 | 5.2 | 10745, 12.5, 13.370, 0.015 |
| V3 raw radii | 0.0151 | 0.0113 | 5.0 | 10745, 12.5, 13.370, 0.015 |
| V6 zero-point map | 0.0110 | 0.0109 | 3.0 | 10546, 11.3, 13.440, 0.035 |

- Martin's model fits worse and shifts the per-run best age by -0.5 Gyr (median; M15: +0.35).
- The zero-point map lowers KS by about 0.004 at every age. Per-run best age shift -0.25 Gyr (median,
  range -1.1 to 0); run 10745 stays at 12.5 Gyr. The V6 age profile is flat to within noise from 11.3 to
  12.5 Gyr with 8 runs, so the full grid (`../tuc_final/`) sets the age.
- Outer colour width (`diag/diag4.py`, upper MS, run 10745): missing spread beyond 1200 px drops from
  about 0.011 to 0.0075-0.008 mag with the map; observed/model width 1.7-1.9 -> 1.4-1.5. The 100-200 px
  bin gets slightly worse (0.0125 -> 0.018). Same residual (~0.008 mag) as M15; still unexplained.

Settings carried to `../tuc_final/setup.json`: binaries 0.009, radial_source raw, zp_correct with
zp_ref 17.9-19.0 and zp_mag_frac 0.08.
