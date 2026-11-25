# tuc_final_mc_indeperr: 47 Tuc with independent per-filter artificial-star errors (2026-11-24)

**Question (the user):** what happens to 47 Tuc with the advisor's error model? In that model the F606W and F814W
artificial-star errors are not taken from the same star and not correlated, as in M. Ying's MakeTestCMD.

**One change** against the Monte Carlo-on fit (`../tuc_final_mc` with production's refinement override, after
`rerun_refine.sh`):

| | this run | `../tuc_final_mc` |
|---|---|---|
| `ast_errors` | "independent": each filter's completeness test and error drawn from its own (radius, magnitude) bins | "nearest": recovery and both errors from one of the 30 nearest artificial stars |

Everything else is identical: isochrones (`tuc_iso_d4v2_mcfix`), the 996 runs, ridge cut, zero-point map, PDMF -0.45,
binaries 0.009, mu 13.15-13.55 and E 0-0.05 (refinement in 0.001), bootstrap at this fit's own lowest cell, and the
statistics.
- `setup.json` is `../tuc_final_mc/setup.json` with only `ast_errors` changed.
- `code/` is identical to `../tuc_final_mc/code`.

**Why the error model matters.** Each row of `47Tucartstars.dat` is one artificial star: one input (x, y) with F606W
and F814W magnitudes, and one output (x, y) with both measured magnitudes. Its two errors correlate (Spearman +0.60 for
the 75,275 recovered stars; a blend brightens both bands). Independent draws lose that correlation. For M15 they
predicted a main sequence 1.7-1.9x wider than observed (`../m15_ast_error_correlation`), and the full M15 run went from
2.2σ to 4.2σ, 14.2 -> 14.6 Gyr (`../m15_final_indeperr`). An 8-run 47 Tuc test on the old isochrones
(`../tuc_errchain`, V0) fitted worse (6.7σ vs 5.2σ) and moved per-run ages by -0.5 Gyr.

**Validation** (`validate_indeperr.py`, `logs/validate_11170062.out`):
1. This code with `../tuc_final_mc/setup.json` reproduces that grid's best cell (10567, 10.2 Gyr, mu 13.50,
   E 0.040): KS 0.0091216 vs 0.0091216.
2. At the joint-error best cells, independent errors give KS 0.0208 (10567) and 0.0227 (10105), against 0.0091.
   That is about 9σ above the null (0.0054 +/- 0.0017), against 2σ. The full grid finds this model's own best cells.

**Driver:** `run.sh` (log `run.log`, `driver.pid`).
1. Wait for `../tuc_final_mc/rerun_refine.sh`.
2. Grid, physics 240 minus my running jobs + 32 free.
3. Refinement.
4. Bootstrap at this fit's own lowest cell.
5. Statistics. `compare_joint_vs_independent.txt` lists the age posteriors of production, joint (Monte Carlo on) and
   independent.

Expected: the refinement rerun ~16:00 on 2026-11-24, this grid ~18 h after, results ~midday 2026-11-25.

## 500 runs (the user, 2026-11-24 15:15)

The user asked for 500 runs of 47 Tuc, done the advisor's way. `runs.txt` is now the first 500 of the 996
(10000-10500); the full list is kept as `runs_996.txt`. The swap was made before `run.sh` reached the grid, so the
driver itself was not edited.
- The grid is about 3,500 units, about 7.5 h after the refinement rerun ends.

`joint500.sh` (log `joint500.log`) does the paired comparison on the same 500 runs:
- After `../tuc_final_mc/rerun_refine.sh`, it takes the joint-error rows of these 500 runs from the corrected
  `../tuc_final_mc_refined/`, checks they have the 0.001 E steps, and runs the statistics in `joint500/` with the joint
  fit's own bootstrap null.
- After `run.sh`, it writes `compare_500.txt`: production (996, Monte Carlo off), joint 996, joint 500 and
  independent 500.

## Grid finished, run 10128 left out (2026-11-25)

**What happened.**
- The grid ran 2026-11-24 16:39 to 2026-11-25 09:18; the cluster was nearly full, so few workers ran at first.
- One grid unit failed: 10128_9.8-10.6 (`logs/11178405_137.out`). For one synthetic sCMD, the ridge traced inside the
  magnitude window came back empty, and `np.interp` failed with "fp and xp are not of the same length"
  (`ridge.py`, `Ridge.__init__`). The same run's unit succeeded with joint errors, so the cause is this draw of
  independent errors, not the isochrone.
- Refinement finished the other 3,499 units. `run.sh` then stopped at `refine.py merge` (`|| exit 1`), because merge
  refuses to run with a unit missing. `joint500.sh` saw the driver end and wrote a comparison without the
  independent rows; that file is kept as `compare_500_joint_only.txt`.

**Fix.** Run 10128 is left out, with no code change:
- `runs.txt` now holds 499 runs; the 500 are in `runs_500.txt`.
- 10128 is a poor fit in both models. Its joint-error best KS is 0.0244, above every joint cut (scaled 0.0175), so it
  carries no weight in `joint500/`. The joint 500-run posteriors are therefore exactly those of these 499 runs, and
  the pairing still holds.
- `resume.sh` (log `resume.log`, PID `resume.pid`) is `run.sh` from the merge on. It also writes `compare_500.txt`,
  with the independent leg labelled `independent_499`.

**Merge.** 22,038 of 30,439 (run, age) cells improved. The best KS is 0.01595 on the grid and 0.01590 refined, at run
10202, 10.6 Gyr, mu 13.44, E 0.004. The joint-error best is 0.00891.

## Results ("all done" 2026-11-25 14:19; `compare_500.txt`, `resume.log`)

**Fit quality.**
- The bootstrap used 10,000 fresh references at this fit's best cell (run 10202, 10.6 Gyr, mu 13.44, E 0.004), drawn
  with independent errors. The null is 0.00516 +/- 0.00182, so the 3σ line is 0.01062. That is the same as the joint
  model's (0.00538 +/- 0.00174): the null depends mostly on the number of stars.
- **The best independent-error fit, KS 0.01590, is 5.9σ above the null mean** (empirical p 0.0002). The joint-error
  best on the same runs is 0.00891, 2.0σ above.
- No row passes 3σ, so strict, strict_D17 and strict_D17conv have no posterior.
- The three best fits disagree on distance and reddening:
  1. run 10202: 10.6 Gyr, mu 13.44, E 0.004;
  2. run 10435: 10.7 Gyr, mu 13.48, E 0.034;
  3. run 10306: 9.8 Gyr, mu 13.55, E 0.047.

**Age posteriors on the same runs (Gyr).** Joint is over the 500 runs. Run 10128 carries no weight there, so this
equals the 499 runs. Each model's shifted and scaled cuts come from its own best fit and null:
- joint: shifted 0.0144, scaled 0.0180;
- independent: shifted 0.0214, scaled 0.0327.

| | joint errors (ours) | independent errors (advisor's) | difference |
|---|---|---|---|
| strict_D17 | 11.84 +/- 0.62 (median 11.69, 8 runs) | none pass 3σ | |
| shifted_D17 | 11.80 +/- 0.69 (11.66, 113) | 11.22 +/- 0.54 (11.06, 19) | -0.58 |
| scaled_D17 | 11.88 +/- 0.77 (11.75, 193) | 11.89 +/- 0.88 (11.76, 243) | +0.01 |
| shifted, whole grid | 11.21 +/- 0.90 (11.09, 147) | 10.71 +/- 0.74 (10.64, 29) | -0.50 |
| scaled, whole grid | 11.26 +/- 0.98 (11.12, 266) | 11.25 +/- 1.10 (11.11, 303) | -0.01 |

**Other parameters, shifted_D17** (joint -> independent):
- [Fe/H] -0.756 -> -0.755;
- mu 13.396 -> 13.410;
- E 0.019 -> 0.016;
- Y 0.256 -> 0.257;
- mixing length 1.70 -> 1.58.

For scaled_D17 every parameter agrees to within 0.02.

**Conclusion.**
- As for M15 (`../m15_final_indeperr`: 2.2σ -> 4-7σ), independent per-filter errors fit 47 Tuc much worse: 5.9σ
  against 2.0σ. Nothing passes the 3σ test.
- The 47 Tuc artificial stars have correlated errors (Spearman +0.60), and drawing each filter separately makes the
  synthetic main sequence too wide.
- The age then depends on how the cut is relaxed:
  - the tight shifted cut keeps 19-29 runs and goes 0.5-0.6 Gyr younger;
  - the loose scaled cut keeps 243-303 runs and does not move.
- Our setup with joint errors stays the better-fitting model.
