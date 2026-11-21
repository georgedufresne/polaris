#!/usr/bin/env python3
"""Did our DSEP4 M15 runs read their Monte Carlo inputs?

DSEP4 (built from ~/dsep4/src) reads the var file (fort.85) through mread only if LMonteCarlo is true
(parmin.f90; default .false. in Common_variables.f90). In the M15 run environment (M15/dsep4/env/nml), only
gs98_3_DATA.nml ([alpha/Fe] = 0) sets LMonteCarlo=.TRUE.; the [alpha/Fe] -0.2/+0.2/+0.4/+0.6/+0.8 control namelists do
not. Runs with [alpha/Fe] != 0 then take KTTAU, FGRY, FGRZ, ALPHAE (no envelope overshoot), ALPHAC and the S-factors from
the physics namelist (phys1.low.nml: KTTAU=5, FGRY=FGRZ=1.0, SStandard central values) and the opacity/neutrino
multipliers stay 1.0. [Fe/H], Y, [alpha/Fe] and the mixing length are passed by the run script, so they vary in all runs.

1. Distribution of each var-file input over our 500 runs, against the value used when it is not read.
2. Ying's isochrones minus ours (same run, 13.4 Gyr, MIST for ours): turnoff colour, SGB-turnoff F606W gap, main-sequence
   and RGB colour; runs that read their var file ([alpha/Fe] = 0, 16 runs) against runs that did not (484).
Output: logs/mc_check_<job>.out, mc_check.csv."""
import sys

import h5py
import numpy as np
import pandas as pd

B = "/dartfs-hpc/scratch/f007058/scmd_tb"
sys.path.insert(0, f"{B}/m15_final/code")
from isochronetoolbox.scmd.model import BCTable, read_isochrones    # noqa: E402

YF = h5py.File("/dartfs-hpc/scratch/f007058/M15/zenodo_ying2025/M15_MC_iso.h5", "r")
mist = BCTable(f"{B}/m15_grid_d4/inputs/bc_mist_acs_wfc_fehm225.npz")
NAMES = ["FeH", "Yprim", "alphafe", "CMIXLA", "FGRY", "FGRZ", "KTTAU", "ALPHAE", "ALPHAC", "S1_pp", "S2_He3He3",
         "S3_He3He4", "S4_pC12", "S5_pC13", "S6_pN14", "S7_pO16", "alexcoef", "opalcoef2", "talphacoef", "plascoef",
         "cocoef", "deboercoef"]
DEFAULT = {"FGRY": 1.0, "FGRZ": 1.0, "KTTAU": 5, "ALPHAE": 0.0, "ALPHAC": 0.0, "S1_pp": 0.9681, "S2_He3He3": 1.0058,
           "S3_He3He4": 1.0, "S4_pC12": 0.8966, "S5_pC13": 1.2727, "S6_pN14": 0.4729, "S7_pO16": 0.9894,
           "alexcoef": 1.0, "opalcoef2": 1.0, "talphacoef": 1.0, "plascoef": 1.0, "cocoef": 1.0}


def frame(g):
    return pd.DataFrame(g["block0_values"][()], columns=[x.decode() for x in g["block0_items"][()]])


def features(v, c):
    up = np.where(v < 5.5)[0]
    i_to = up[np.argmin(c[up])]
    post = np.where((np.arange(len(c)) > i_to) & (c >= c[i_to] + 0.05))[0]
    i_sgb = post[0]
    ms = np.where(np.arange(len(c)) < i_to)[0]
    rgb = np.where(np.arange(len(c)) > i_sgb)[0]
    vs = v[i_sgb]
    return dict(c_to=c[i_to], gap=vs - v[i_to], c_ms45=np.interp(4.5, v[ms][::-1], c[ms][::-1]),
                c_rgb=np.interp(vs - 1.5, v[rgb][::-1], c[rgb][::-1]) if len(rgb) and v[rgb].min() < vs - 1.5 else np.nan)


var = pd.DataFrame([[float(l.split()[0]) for l in open(f"/dartfs-hpc/scratch/f007058/M15/dsep4/var/varfeh-227.{r}")]
                    for r in range(10000, 10500)], columns=NAMES, index=range(10000, 10500))
read = pd.Series(np.isclose(var.alphafe, 0.0), index=var.index)
print(f"runs that read their var file ([alpha/Fe] = 0): {read.sum()} of {len(var)}")
print("\n1. var-file inputs over our 500 runs (median, 16-84%) vs the value used when not read:")
for k, d in DEFAULT.items():
    x = var[k]
    extra = f"; values {dict(x.value_counts().sort_index())}" if k == "KTTAU" else ""
    print(f"   {k:11s} median {x.median():.4f} ({x.quantile(.16):.4f}..{x.quantile(.84):.4f}); not read -> {d}{extra}")

rows = []
for r in var.index:
    iso = read_isochrones(f"/dartfs-hpc/scratch/f007058/M15/iso_dsep4/{r}/isos/isochrones_example.txt")
    g = iso[np.isclose(iso.isochrone_age_yr, 13.4e9)].sort_values("EEP")
    if not len(g):
        continue
    m = mist.absolute(g)
    try:
        fo = features(m["F606W"], m["F606W"] - m["F814W"])
        d = frame(YF[f"{r}/isochrones/cmd/13400"]).sort_values("EEP")
        fy = features(d.F606W.values, d.F814W.values)
    except (IndexError, KeyError):
        continue
    rows.append(dict(simulation_id=r, read=bool(read[r]), KTTAU=var.KTTAU[r], ALPHAC_ying=frame(YF[f"{r}/var"]).ALPHAC[0],
                     **{f"d_{k}": fy[k] - fo[k] for k in fy}))
t = pd.DataFrame(rows)
t.to_csv("/dartfs-hpc/scratch/f007058/scmd_tb/m15_dsep4_montecarlo_check/mc_check.csv", index=False)
print(f"\n2. Ying minus ours at 13.4 Gyr (ours with MIST), median (16-84%) [n]:")
for name, s in (("read var file ([a/Fe]=0)", t[t.read]), ("did not read ([a/Fe]!=0)", t[~t.read])):
    print(f"   {name:26s} " + "; ".join(f"{k} {s[f'd_{k}'].median():+.4f} ({s[f'd_{k}'].quantile(.16):+.4f}..{s[f'd_{k}'].quantile(.84):+.4f})"
                                            for k in ("c_to", "gap", "c_ms45", "c_rgb")) + f" [{len(s)}]")
print("\n   runs that did not read, by the KTTAU their var file asked for:")
for k, s in t[~t.read].groupby("KTTAU"):
    print(f"   KTTAU {int(k)}: gap {s.d_gap.median():+.4f}, turnoff colour {s.d_c_to.median():+.4f}, MS colour {s.d_c_ms45.median():+.4f} [{len(s)}]")
s = t[t.read]
print(f"   runs that read: Spearman of gap difference with his ALPHAC {s.d_gap.corr(s.ALPHAC_ying, method='spearman'):+.2f}")
