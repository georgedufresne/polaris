"""Summarise the orientation audit: data vs model and split-AST null per ECDF orientation."""
import glob
import json
import sys

import numpy as np
import pandas as pd

H = "/dartfs-hpc/scratch/f007058/scmd_tb/bootstrap_audit/orientation"
SPLIT = "/dartfs-hpc/scratch/f007058/scmd_tb/tuc_grid_nearest_bin02_split/bootstrap_2d_ks_results.csv"


def load(pattern):
    rows = []
    for p in sorted(glob.glob(f"{H}/{pattern}")):
        rows += [json.loads(line) for line in open(p) if line.strip()]
    return pd.DataFrame(rows)


data, null, grid = load("data_*.jsonl"), load("null_*.jsonl"), load("grid.jsonl")
print(f"data realisations {len(data)} (ids {sorted(data.id.tolist())}), null realisations {len(null)}, grid rows {len(grid)}")

# ---- validation ------------------------------------------------------------------------------
sp = pd.read_csv(SPLIT).set_index("bootstrap_number")
if len(null):
    m = null.set_index("id")["bb_p95"]
    diff = (m - sp.loc[m.index, "ks_statistic"]).abs()
    print(f"null bb_p95 vs split-null CSV for the same {len(m)} ids: max |diff| {diff.max():.3e}")
    sub = sp.loc[m.index, "ks_statistic"]
    print(f"  CSV, these ids: {sub.mean():.5f} +- {sub.std():.5f};  CSV, all 10,000: "
          f"{sp.ks_statistic.mean():.5f} +- {sp.ks_statistic.std(ddof=0):.5f}")
for df, name in ((data, "data"), (grid, "grid"), (null, "null")):
    if "ks_pipeline" in df:
        v = df.dropna(subset=["ks_pipeline"])
        print(f"{name}: pipeline ks_2d vs bb_p95 max |diff| {np.abs(v.ks_pipeline - v.bb_p95).max():.3e}")
    for k in ("bb", "bf", "rb", "rf"):
        c = f"{k}_ex_identity_maxdiff"
        if c in df:
            print(f"{name}: exact-ECDF identity vs literal masked rank2d, {k}: max |count diff| {df[c].max()}")
for df, name in ((data, "data"), (grid, "grid"), (null, "null")):
    if len(df):
        print(f"{name}: cross ties colour {df.cross_ties_c.max()}, magnitude {df.cross_ties_m.max()}; "
              f"obs points with interpolated sCMD ECDF = 0: " +
              ", ".join(f"{k} {int(df[f'{k}_nfi0'].mean())}" for k in ("bb", "bf", "rb", "rf")))
if len(grid):
    print(f"grid sCMD: pipeline ks_2d {grid.ks_pipeline.iloc[0]:.10f} (grid CSV 0.0168281775635039)")

# ---- table -----------------------------------------------------------------------------------
rng = np.random.default_rng(1)


def sep_err(d, n, nboot=4000):
    """Bootstrap (over realisations) standard error of (mean_d - mean_n) / std_n."""
    d, n = np.asarray(d), np.asarray(n)
    out = []
    for _ in range(nboot):
        dd = d[rng.integers(0, len(d), len(d))]
        nn = n[rng.integers(0, len(n), len(n))]
        out.append((dd.mean() - nn.mean()) / nn.std(ddof=1))
    return float(np.std(out))


rows = []
stats = []
for base in ("p95", "p99", "ex_p95", "ex_p99", "z_p95", "z_p99", "ex_max"):
    for k in ("bb", "bf", "rb", "rf"):
        stats.append((f"{k}_{base}", k, base))
    stats.append((f"ffmax_{base}", "max4(stat)", base))
for base in ("p95", "p99", "ex_p95", "ex_p99", "ex_max"):
    stats.append((f"ffpt_{base}", "pointmax4", base))
stats += [("ying_p95", "Ying(bf,renorm,mask0)", "p95"), ("ying_p99", "Ying(bf,renorm,mask0)", "p99"),
          ("c1d_p95", "1D colour", "p95"), ("m1d_p95", "1D F606W", "p95"),
          ("c1d_max", "1D colour", "max"), ("m1d_max", "1D F606W", "max")]
for col, orient, base in stats:
    if col not in data or col not in null:
        continue
    d, n = data[col].values, null[col].values
    nm, ns = n.mean(), n.std(ddof=1)
    lim = nm + 3 * ns
    g = float(grid[col].iloc[0]) if len(grid) and col in grid else np.nan
    rows.append(dict(stat=base, orientation=orient, data_mean=d.mean(), data_std=d.std(ddof=1),
                     null_mean=nm, null_std=ns, null_3sig=lim, sep=(d.mean() - nm) / ns,
                     sep_se=sep_err(d, n), ratio=d.mean() / nm, data_pass=f"{int((d <= lim).sum())}/{len(d)}",
                     grid=g, grid_sep=(g - nm) / ns, null_max=n.max()))
T = pd.DataFrame(rows)
pd.set_option("display.width", 250)
pd.set_option("display.max_rows", 200)
fmt = {c: "{:.5f}".format for c in ("data_mean", "data_std", "null_mean", "null_std", "null_3sig", "grid", "null_max")}
fmt.update(sep="{:.1f}".format, sep_se="{:.1f}".format, grid_sep="{:.1f}".format, ratio="{:.2f}".format)
print(T.to_string(index=False, formatters=fmt))
T.to_csv(f"{H}/summary_table.csv", index=False)
print(f"\nwritten {H}/summary_table.csv")
