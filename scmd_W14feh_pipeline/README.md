# NGC 6791 sCMD leg: W v1.4 with [Fe/H] converted to our Z convention ("W14feh", 2026-11-22)

## Why

The W training grid computed Z with Martin's `gs98getz`, which gives exactly 1.121 x ours. Every one of our DSEP grids
(OG, DSEP3 physics-on "New", and the DSEP4 grids) uses ours. `Isochrones_W14` passed each var file's [Fe/H] straight to
W, so W modelled 12 % more Z than the DSEP3 run of the same var file. As a result its turnoff was 0.002 dex cooler
(`../emul_W14/feh_convention`).

The user asked (2026-11-22) to rerun the whole W v1.4 leg with the matching convention and re-age every leg.

## What changes

**Only the isochrones.** `../Isochrones_W14feh/<run>/` comes from `../emul_W14/gen_W14_iso.py --z-factor 1.121`
(batch `../emul_W14/gen_W14feh.batch`):
- W and its age head get [Fe/H]_W = getz^-1(getz([Fe/H]) / 1.121), about -0.05 dex;
- the header keeps the var [Fe/H] and Z = getz([Fe/H]).
- Checked on run 20041: the W [Fe/H] matches the convention test, and the turnoff moved to within 0.0004-0.0007 dex of
  DSEP3 (from +0.002).

The run list is the 668 KTTAU 0/5 runs of the W14 leg (`../emul_W14/runs_W14feh.txt`).

**The scripts here** are copies of `../scmd_W14_pipeline`, changed only in paths:
- the new folders: `Isochrones_W14feh`, `sCMD_W14feh`, `2d_ks_W14feh`, `plots_for_2d_W14feh`, `OUTPUT_KDE_W14feh`,
  `analysis_W14feh`;
- the home inputs, now read from md5-checked copies in `inputs/`.

**The bolometric-correction table.** The original picks `os.listdir(~/.fidanka/bol/MIST/UBVRIplus)[3]`, which in home
is `fehm300.UBVRIplus`, the [Fe/H] = -3.0 table, for every NGC 6791 leg. That folder has not changed since 2025-07.
The copy loads that same file from `inputs/`. The interpolated G, BP and RP corrections are identical (max difference
0.0 at 1,200 points).

This is a pre-existing quirk shared by all legs, so the comparison is unaffected. It does matter for the absolute
NGC 6791 ages.

**Python and home access.** Python comes from the scratch env, and `HOME` is set to scratch, so no batch job touches the
Kerberos home mount. Generation is the exception: it needs torch, so it uses the home conda env (short jobs).

**Smoke test.** One age of run 20041 went through sCMD + KDE; the outputs landed in the new folders and were pruned to
4 files.

## Driver: `run_W14feh.sh` (log `run_W14feh.log`)

1. Isochrones, with one retry.
2. sCMD array, throttled to 240 physics CPUs minus my other jobs, with finished-but-hung tasks reaped.
3. Wait for the DSEP3 (11131955) and W14 (11131934) reruns.
4. Completeness: every KTTAU 0/5 run in the New, W14 and W14feh legs needs all of its isochrone's ages (81) in its
   2D-KS table. Incomplete tables are moved to `<ksdir>_partial/` and rerun once.
5. Merge + prefilter for all three legs. The previous `filtered_ks_scmd.csv` is kept as
   `filtered_ks_scmd_pre_rerun.csv`.
6. `scmd_6way_compare.py` (-> `../scmd_6way_comparison.txt`), `age_breakdown_W14feh.py` and the W14 `age_breakdown.py`.

Then `../scmd_tb/tuc_final_mc/run.sh` starts the 47 Tuc refit.

## Result (driver "all done" 2026-11-23 10:13)

**Completeness.** W14feh 668/668, New 0 incomplete. The W14 leg is 667: run 20374 has 80 ages, because its
`Isochrones_W14` 7.2 Gyr curve is truncated (`../emul_W14/iso_truncation`). It was retried once and kept out.

**Posterior ages (KTTAU 0/5, each leg's own sets)** (`../scmd_6way_comparison.txt`):

| leg | age (Gyr) |
|---|---|
| OG | 8.496 +/- 0.620 |
| New | 8.534 +/- 0.647 |
| W14 | 8.363 +/- 0.647 |
| W14feh | 8.390 +/- 0.668 |

**Paired difference on common sets (bootstrap SE):**

| pair | difference (Gyr) | n |
|---|---|---|
| W14feh - New | -0.139 +/- 0.017 | 374 |
| W14 - New | -0.188 +/- 0.019 | 354 |
| W14feh - OG | -0.075 +/- 0.028 | |
| W14feh - W14 | +0.049 +/- 0.022 | |
| New - OG | +0.051 +/- 0.027 | |

**Decomposition of the pooled gap** (`age_gap_decomposition.py` -> `age_gap_decomposition.txt`):

| | pooled | within-set shift | set reweighting |
|---|---|---|---|
| W14 - New | -0.188 | -0.239 | +0.051 |
| W14feh - New | -0.139 | -0.078 | -0.061 |

- The Z fix removes about two-thirds of the per-set age offset. The per-set median moved from -0.24 to -0.05 Gyr
  (`age_breakdown.txt`), and the He correlation of the per-set gap dropped from rho +0.29 to +0.15.
- The Z-matched W gives relatively more posterior weight (more and better-KS rows) to the sets that DSEP3 puts young
  (corr -0.25). That is why the pooled gap shrinks by only 0.05 Gyr.
- Fit quality is equal: median lowest KS 0.0306 (W14feh) vs 0.0303 (New).

**All parameters, mean vs median, and what is left of the gap** (`param_gap_breakdown.py` -> `param_gap_breakdown.txt`,
2026-11-23):

- **Age posteriors, own sets (mean / median, Gyr):** OG 8.496 / 8.466, New 8.534 / 8.498, W14 8.363 / 8.319, W14feh
  8.390 / 8.339. Medians sit 0.03-0.05 below the means.
- **W14feh - New, age differences (Gyr):**

  | measure | mean | median |
  |---|---|---|
  | posterior | -0.139 +/- 0.017 | -0.155 +/- 0.023 |
  | per set | -0.064 +/- 0.016 | -0.049 +/- 0.016 |

- **Other parameters, W14feh - New:** DM +0.009, E(B-V) -0.002, [Fe/H] +0.008, Y -0.0004, [alpha/Fe] +0.008, mixing
  length -0.014. Each is at most 0.2 of the posterior sd.
- **Within-set age shift (-0.078 +/- 0.013).** In each set, W lands slightly farther (+0.008 mag) and redder (+0.0025).
  The per-set age shift at equal DM and E(B-V) is +0.03 Gyr, so the age goes with the distance/reddening degeneracy.
  `../emul_W14/feh_convention` (W fitted to DSEP3's own isochrones, distance + reddening free) gives the same -0.05 to
  -0.10 Gyr. It is an emulator-accuracy effect, not a pipeline one.
- **Reweighting (-0.061 +/- 0.013).** All of it goes through the var-file inputs (residual +0.004). W14feh favours sets
  with slightly higher [Fe/H] (-0.020 Gyr), lower mixing length (-0.020), higher [alpha/Fe] (-0.015) and alexcoef
  (-0.009), which DSEP3 ages young. This term was +0.051 before the Z fix.
