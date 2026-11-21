# The DSEP4 side of these figures did not read its Monte Carlo inputs (found 2026-11-19)

`plot_d3_d4_isochrones.py`'s docstring says the DSEP4 isochrones (`../tuc_iso_d4v2`) had "Monte Carlo ON". They did not:
- the 47 Tuc DSEP4 tracks (`tuc_dsep4_batch/run_tuc_dsep4_mc.batch`, `NML=~/nml_dsep4_mc`) used namelists that set
  `LMonteCarlo` only for [alpha/Fe] = 0;
- all 1000 var sets (`~/mcdsep/MC_var_tuc_d4`) have [alpha/Fe] 0.2 (399) or 0.4 (601);
- DSEP4 defaults `LMonteCarlo` to false, so no run called `mread`.

The logs confirm it. The message "triple-alpha luminosity reached 1% of total", printed only inside `if (lmontecarlo)`, is
absent from the 47 Tuc DSEP4 logs (runs 10000, 10005, 10010, 10091). It is present for an M15 [alpha/Fe] = 0 run and for
the fixed M15 rerun.

DSEP3 (`../tuc_iso`, `~/dsepX_fixed`) did read its var files (`bd_montec.f`: `DATA LMonteCarlo/.TRUE./`). So in
`d3_d4_*.png`, the DSEP4-DSEP3 differences combine the code version with DSEP4's missing Monte Carlo physics.

**Redrawn** (`plot_20_panels_mcfix.py` -> `mcfix/`): DSEP4 from `../tuc_iso_d4v2_mcfix`. These are the figure's 20 runs
rerun with `LMonteCarlo` on (`tuc_dsep4_batch/run_tuc_dsep4_mcfix.batch`, `nml_dsep4_mc_fixed`) and rebuilt with the same
v2 builder. `_with_old` versions also show the old physics-off DSEP4 in grey.

Also note: the 47 Tuc production fit (`../tuc_final`) uses `../tuc_iso_d4v2`.

## Result (redrawn 2026-11-20 09:21, `mcfix/d3_d4_20_panels_{hrd,cmd}[_with_old].png`)

With both models reading their Monte Carlo inputs, DSEP4 lies on top of DSEP3 in essentially every panel. Shown as
DSEP4 - DSEP3 in log Teff (1e-3 dex), median (90% of |values| below), over the 20 runs at 8 and 12 Gyr:

| where | DSEP4 physics on | old DSEP4 (physics off) |
|---|---|---|
| lower MS (log L -0.5) | -0.4 (1.9) | +1.0 (5.4) |
| upper MS (log L 0.0) | +0.1 (1.4) | +1.1 (5.6) |
| turnoff | +0.2 (1.6) | +1.5 (5.9) |
| RGB (log L 1.0) | +0.0 (1.1) | +5.4 (14.6) |
| RGB (log L 1.5) | +0.1 (1.2) | +4.6 (14.9) |

- **The code difference is negligible.** DSEP3 and DSEP4 with the physics on agree to about 1e-3 dex (~15 K) everywhere.
- **The old figure's differences came from the missing Monte Carlo physics.** They were 5-10x larger, mostly a hotter
  RGB, from the default surface boundary condition and the missing diffusion and overshoot variations.
