# M15 with 4,000,000-star model sCMDs (2026-11-11)

Purpose: Ying et al. (2025, section 4) build 4,000,000 points per isochrone. Every production sCMD has
exactly as many stars as the observed fit sample (32,085), so the model's ECDF is about as noisy as the
data's. This experiment changes only the model size, for two baselines:

| run | baseline | zero-point map | E range |
|---|---|---|---|
| `prod/` | `../m15_final_refined/` (production) | on | 0.000-0.120 |
| `nozp/` | `../m15_yingbox_nozp_refined/` | off | 0.080-0.150 (Ying's box, applied directly) |

Both use mu 15.30-15.60. Each `setup.json` is its baseline's refined setup plus `n_model` 4,000,000 and
`ks_model_ecdf` "count".

## What changes, and what does not

- **Model size:** `n_model` (new in `code/`) sets the model sCMD size. The ridge of each sCMD is still traced
  through an observed-sized sample (as many stars as the observed trace box), as in production. The sCMD is
  then filled to 4M stars inside that ridge band, so the selection is unchanged.
- **Model ECDF:** with millions of points it is counted exactly at each observed star (`stats.ecdf_at`)
  instead of linearly interpolated. That is the same quantity without the interpolation error.
- **Bootstrap:** fresh references keep the observed number of stars and use all artificial stars; each
  realisation's model has 4M stars and resampled artificial stars. This is data-sized against 4M, as in the
  fit (and as Ying's null, data-sized against 1M).
- **Unchanged:** everything else. That covers the crn replay, the 95th-percentile KS, the nearest-AST
  errors, the radii, and the strict / shifted / scaled posteriors.

## Search

The production grid cannot be rerun at 4M stars per cell (about 18 s each). `refine4m.py` therefore
descends from the 32k refined minimum of every (run, age) whose 32k KS is at most 2.5x the grid's best. It
uses the 8-neighbour descent of `ks_refine/refine.py` on the same lattice (mu step 0.01, E step 0.001).

| run | (run, age) starts | runs | 32k KS limit |
|---|---|---|---|
| prod | 11,334 | 497 | <= 0.0342 |
| nozp | 10,627 | 426 | <= 0.0474 |

The scaled posterior's cut is roughly 2x the best KS, so the 2.5x limit should cover it.

## Files

- **`code/`:** copy of `../m15_final/code`. The changes are `Setup.n_model` (0 reproduces production),
  `Setup.ks_model_ecdf` ("interp" reproduces production), `stats.ecdf_at` / `ks_2d_count`, and the
  bootstrap reference at the observed size.
- **Checks:**
  - `time_cell.py`: timing at 32k, 250k and 1M stars.
  - `validate_count.py` (`logs/validate_count.out`): the exact count against brute force and against the
    interpolated ECDF.
  - `test/`: the descent replays the production best cell with 32k sCMDs and stops there.
  - `test4m/`: smoke test of one 4M evaluation.
- **Search:** `refine4m.py`, `refine4m.sbatch`, `prod/starts.csv`, `nozp/starts.csv`.
- **Driver:** `run.sh` (`driver.pid`, log `run.log`). It waits for `../tuc_yingbox_nozp/run.sh`, then runs
  the descent (nozp first), merge, bootstraps, statistics and figures.
- **Outputs:**
  - `prod/` and `nozp/`: `combined_2d_ks_results.csv` (4M cells only), `bootstrap_2d_ks_results.csv`,
    `posterior_*`;
  - `analysis/ks_vs_age_32k_vs_4m.png`.

## Validation

1. **Exact count:** `ecdf_at` equals a brute-force count on random sets, ties included (difference 0).
2. **Interpolated vs counted KS** on production's best cell (run 10357, 14.2 Gyr, mu 15.58, E 0.096). The
   difference shrinks as the model grows:

| model sCMD | KS interpolated | KS counted | count time vs interpolation time |
|---|---|---|---|
| 32k | 0.013690 | 0.013869 | 0.2 s vs 1.8 s |
| 250k | 0.014352 | 0.014466 | 0.3 s vs 12.1 s |
| 1M | 0.014480 | 0.014549 | 0.6 s vs 48.1 s |

3. **Replay (32k, interpolated):** the descent starts at production's best cell with KS 0.0136905, as in the
   file. All 8 neighbours are higher, so it stops. The neighbours one step away (0.01 in mu, 0.001 in E)
   have KS 0.0141-0.0196: with 32k models the KS surface is jagged at the lattice scale.
4. **4M cells:** at production's three reference cells, each takes 15 s to generate and 3 s to score, and
   peaks at 1.2 GB. Production's best cell rises from 0.0137 (32k) to 0.0148 (4M). The 32k minimum was
   partly a favourable noise draw.

## Model-size check at the age-deciding cells (`model_size_check.py`, `logs/model_size_check.out`)

Run on 2026-11-11, before the 4M run. Model sCMDs of 1x, 2x, 4x and 16x the observed count, with 3-5 seeds
each and the counted ECDF at every size.

| cell | 1x | 2x | 4x | 16x |
|---|---|---|---|---|
| production, best (14.2 Gyr) | 0.0158 +/- 0.0011 | 0.0157 +/- 0.0006 | 0.0154 +/- 0.0004 | 0.0152 +/- 0.0002 |
| production, 12.3 Gyr | 0.0228 +/- 0.0012 | 0.0223 +/- 0.0007 | 0.0219 +/- 0.0003 | 0.0224 +/- 0.0003 |
| no map, best (14.0 Gyr) | 0.0206 +/- 0.0007 | 0.0202 +/- 0.0004 | 0.0200 +/- 0.0001 | 0.0200 +/- 0.0001 |
| no map, 12.3 Gyr | 0.0209 +/- 0.0006 | 0.0204 +/- 0.0003 | 0.0203 +/- 0.0003 | 0.0208 +/- 0.0001 |

The +/- figures are the std across seeds. KS(12.3 Gyr) - KS(best), mean over seeds:

| | 1x | 2x | 4x | 16x |
|---|---|---|---|---|
| production | +0.0070 | +0.0067 | +0.0065 | +0.0072 |
| no map | +0.0002 | +0.0002 | +0.0004 | +0.0007 |

- **Bigger models cut the noise, not the KS level.** The seed-to-seed std drops by 4-10x from 1x to 16x,
  but the mean KS drops by only 0.0001-0.0006.
- **Production's grid value was a lucky draw.** Its 0.0137 (seed 0 with the interpolated ECDF; 0.0139
  counted) sits about 1.7 sigma below the 1x mean of 0.0158. The grid minimum is biased low by selection.
- **Which age wins does not change.** Production prefers 14.2 over 12.3 Gyr by 0.007 at every size. Without
  the map the preference for 14.0 over 12.3 grows, from 0.0002 (within the noise) to 0.0007 (7x its noise).
  The bigger model does not push toward a younger age.
- **Duplicating the observed stars changes nothing:** KS 0.013869 either way. Every star counted twice
  leaves the ECDF the same.

These are fixed cells from two runs, not re-optimised. The 4M descent covers every run.

## Results

Not run. The user cancelled the 4M descent on 2026-11-11, while the driver was still waiting for 47 Tuc and
before any job was submitted (the driver was stopped; `run.log` records it). The model-size check above is the
outcome: bigger models cut the seed noise but do not change which age wins.
