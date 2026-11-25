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
