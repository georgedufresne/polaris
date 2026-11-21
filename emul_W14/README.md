# NGC 6791 isochrones from W v1.4 (2026-11-18)

Purpose: the user updated the W emulator (CONF1DENCE-W, github.com/georgedufresne/CONF1DENCE-W, commit fcfac3c
"Add W v1.4 (DSEE v2, pyeep-1.4 EEPs) and make it the default"). This regenerates the NGC 6791 W isochrones with it and
runs them through the same sCMD pipeline as the October W leg (`../emul_W`, `../Isochrones_W`, `../scmd_emul_pipeline`),
to compare with the other NGC 6791 sCMD legs (OG, New, Emul(old), W; `../scmd_4way_comparison.txt`).

## Inputs (copies on scratch, so batch jobs do not need the Kerberos home mount)

- `CONF1DENCE-W_fcfac3c/`: snapshot of the home clone after `git merge --ff-only origin/main` (ec8b4d9 -> fcfac3c; there
  were no local changes). `COMMIT` records the hash. The W v1.4 weights have the same md5 as the home copy.
- `var/`: the 1000 NGC 6791 var files `varfeh+031.20000-20999` (from `~/mcdsep/MC_var_noco`, as the October run).
- `gs98getz`: from `~/mcdsep/run`.
- Python: the conda env `fidanka-env` (`~/.conda/envs/fidanka-env`: torch 2.8.0, zuko 1.2.0, as W requires), the same
  env as the October run.

## gen_W14_iso.py

This is a copy of `../emul_W/gen_W_iso.py` with the v1.4 recipe of `isochrone_W.py`:
- pyeep-1.4 EEP grid (202-1055 = ZAMS to RGB tip, 2,695 EEPs);
- `tol = 5`, with tip-anchored monotonic mass;
- the flow's mode for the positions;
- the age table built once per run (it does not depend on the age).

Kept from October: ages 4.0-12.0 Gyr in 0.1, masses 0.10-2.5 Msun on 1,200 points, the var-file -> label mapping
(C12+He4 = 1.0) and the MIST-format output to `../Isochrones_W14/<run>/isos/isochrones_example.txt`.

Conductive opacity: W v1.4 takes a perturbation in [-1, 1] (default 0; the DSEP4 Blouin20 convention,
`dsep4/src/mainopac.f90`). These varfeh files hold a multiplier (cocoef ~1 +/- 0.2). `--cond var` (the default) passes
the var value as October did; `--cond zero` uses 0. The validation measures the difference. The NGC 6791 2D KS uses a
1.7 mag box around the turnoff (`../scmd_why_no_change.txt`), where conductive opacity acts little.

## Validation (`validate.batch`, `logs/validate_<job>.out`)

`gen_W14_iso.py --validate 20000` compares this generator with `isochrone_W.w_isochrone` at 4.0, 8.0, 8.3 and 12.0 Gyr
(they must be identical), measures the conductive-opacity effect, and times one full run (81 ages, `test_out/`).

## NGC 6791 DSEP grids: no rerun needed

The M15 DSEP4 grid had a Monte Carlo read bug (`../scmd_tb/m15_dsep4_montecarlo_check`). The NGC 6791 "New" (noco)
grid was built with `~/dsepX_fixed` (`~/mcdsep/run/runmc2_noco.sh`). There `bd_montec.f` sets
`DATA LMonteCarlo/.TRUE./`, so `mread` reads every var file. "OG" is the older ChaboyerB grid that "New" supersedes
(`../scmd_why_no_change.txt`).
