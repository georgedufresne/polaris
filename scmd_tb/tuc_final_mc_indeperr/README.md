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
