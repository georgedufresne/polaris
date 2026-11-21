## Bugs found

1. **NGC 6791's `ks_3sigma = 0.0468` carries over into the 47 Tuc setup** (pipeline.py:70). As a result `filtered_ks_scmd.csv` keeps 3,464 of 23,668 rows, but none of them pass either 47 Tuc null. The posterior is not affected, because it uses the bootstrap's own limit.
2. **NGC 6791's KDE bandwidth 0.15 is hard-coded** (stats.py:31, :97). It chooses the 4 cells per (run, age) that go on to the KS stage, while the 47 Tuc MS is only 0.02–0.03 wide. The KS-best cell is known to make it through only at run 10099, 11.4 Gyr, so better cells may never have been KS-scored (untested).
3. **Wrong docstrings:** `rank2d` does not count strictly when there are ties (KS effect about 0.0002). The calibration notes in bootstrap.py assume that comparing half A against half B is the right test.
4. **Fixed split:** the AST split always uses `RandomState(0)`, and that split gives one of the wider nulls of the splits tried.

There is no numerical bug in the bootstrap code itself.

## 1. How the bootstrap works

1. **Cell:** the one given by `--cell`, otherwise the grid's lowest-KS row (run 10099, 11.4 Gyr, μ 13.38, E(V−I) 0.015, KS 0.01683).
2. **Split:** `RandomState(0)` divides the 99,999 ASTs into two halves with no overlap: A (50,054) and B (49,945). Each half is rebuilt with its own neighbour pools and centre. `--shared-ast` uses the full list everywhere instead.
3. **Reference:** using half A, the grid's random stream is replayed through 257 cells. That stream does not depend on the AST list or the ridge. So the reference reuses the random numbers of the grid's sCMD (masses, binaries, radii, neighbour picks), but gets half-A errors, losses and ridge cut. It is a single fixed draw, and its ECDF is rank2d/N.
4. **Realisations** i < 10,000: using half B, `RandomState(i)` makes a fresh sCMD of the same cell with the grid's generator and its own ridge cut, keeping 47,413 stars.
5. **Score:** the same `ks_2d` used for the data, with the reference standing in for the data. The realisation's ECDF is interpolated at the reference points, and KS = 95th percentile of |ΔECDF|.
6. **Collect:** split 0.008873 ± 0.001920, giving mean + 3σ = 0.014633; shared 0.003948 ± 0.001162, giving 0.007435. The posterior applies this one cell's limit to every row.

**Shared between reference and realisations:** cell, generator, observed radii, ridge settings, N, the statistic, and the split (with `--shared-ast`, also the AST list). **Independent:** the star draws, the ridge traces, and in split mode the AST lists.

## 2. What differs between data-vs-model and the bootstrap

| Asymmetry | Direction | Size |
|---|---|---|
| AST-list noise: in data vs model, only the model copies errors from a finite list (the full one) | shared null: none (too narrow); split: two half-lists, nominally 4× the variance (too wide) | full/full 0.0043, half/full 0.0053 (the matched case), half/half 0.0080 |
| One reference, one split | null depends on both draws | seed-0 reference means 0.0068–0.0121 (3σ 0.0103–0.0175); splits 1 and 2 give ≈0.0055 |
| Best cell = minimum of 23,668 rows | lowers data KS | 0.0168 vs 0.01846 ± 0.00094 (48 fresh realisations) |
| Grid node vs true parameters | raises data KS | ≤0.0005 |
| Interpolated model ECDF | lowers data KS | 0.0005–0.0007 (null ≤0.00014) |
| Ties, hull, star counts, ridge, code version | negligible | ≤0.0002 |
| Model radii taken from already-detected stars (losses applied twice); few ASTs in the core | probably raise data KS only | core under-weighted by ~10–15%; KS effect not measured |
| Missing physics (colour offsets, blue wing, wide core MS, multiple populations, 2% binaries, field stars) | raises data KS only | residual ≈0.017; which physics is responsible is untested |

## 3. What we tested

| Hypothesis | Test | Key numbers | Verdict | Verification |
|---|---|---|---|---|
| Counting bluer-and-brighter stars inflates the data KS | 4 quadrants, Fasano–Franceschini, Ying's recipe; 12 model vs 40 null realisations | vs split null: pipeline 5.4σ, Ying's quadrant 9.4σ (his full recipe 6.3σ), other quadrants 13.8–25.9σ; 0/12 pass; colour and V marginals alone 4.3σ, 7.2σ | does not explain (the pipeline's quadrant is the least sensitive) | not verified |
| Ying-style noise-free model ECDF changes the outcome | 1,137,912-star pooled model, exact ECDF | data 0.01779, null 0.00771 ± 0.00156 (3σ 0.01233): 6.5σ, vs 6.3σ comparing equal-size samples | does not explain | not verified |
| Split null counts AST noise twice over (4× vs 1×) | code check; combining the two tables; 835 KS values from mirrored designs | mirrored null 0.00519 ± 0.00174, 3σ 0.0104; best fit 6.7σ (4.1σ against the split null) | explains a large part of why the fit looked close; fixing it makes the test stricter | reproduced twice (3σ 0.0097, z 8.2; non-overlapping lists 3σ 0.0116); the 4×-vs-1× scaling is only approximate |
| Null misses the grid-node penalty | fits to grid tables; Monte Carlo | ≤0.0005 at the best node; split 3σ with the penalty 0.0150 | small part | not verified |
| A single reference and winner's curse distort the comparison | 10 references × 48 realisations; 48 fresh data realisations; empirical Bayes | reference sd 0.00166; the official reference is typical; grid 0.0168 is a low draw at 1.73σ; corrected ≈0.0175 | does not explain; both effects favour the model | not verified |

## 4. Bottom line

- **Nothing in the bootstrap or KS machinery explains the KS values.** Its effects (section 2) are small, and taken together they favour the model.
- **The best fit only looks close to passing** because the split null counts AST-list noise several times over, seed 0 happens to be a wide split, and the grid minimum is a lucky low draw.
- **The best-justified null** takes reference errors from one AST half and realisation errors from the full list, pooled over references and splits. It has mean ≈0.0052–0.0056 and 3σ ≈0.0097–0.0116.
- **Where the best models sit:** the raw 0.0168 is about 6–9σ above that null, and the corrected ≈0.0175 is further still.
  - The mismatch (≈0.017, i.e. 1.7% in the ECDF) is small, but with 47,413 stars the test detects it easily.
  - It sits mostly at V 17.3–18.0 (turnoff and upper MS), where the data have ≈2% more stars at and above the turnoff.
- **To pass the best-justified null, KS must drop to ≈0.0097–0.0116.** Even against the split limit, the corrected value would have to drop by 0.0021–0.0039.
  - *Justified corrections:* the mirrored null, pooling over references and splits, the corrected KS, exact ECDFs, a discretisation term, and scoring the cells stage 2 skipped. Together they push the fit *further from* passing. Fits predict 0.0152–0.0168 for the skipped cells (untested), which would still fail.
  - *Changes that only loosen the test:* keeping the split null; picking a lenient reference (R9 accepts 0.0168 but not 0.0175) or split; inflating the null std 2.2× (to 0.0034); fitting only ~10,000–20,000 stars (a trace estimate, untested). No quadrant, percentile or statistic we tried makes it pass.
  - Only changes to the model can honestly close the gap. A local colour correction is reported to give 0.0140 (not re-tested here); that passes 0.0146 but not the mirrored limit.
- **Untested:** other cells and runs; using one cell's null for all rows; the KS effect of the core under-weighting, of the sparse core ASTs, and of the missing physics.

All test scripts and results are in /dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/