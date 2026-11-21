#!/usr/bin/env python3
"""Run 20000: dsepX_fixed vs dsep4 with Monte Carlo mode off (as first run) vs dsep4
with Monte Carlo mode on (LMonteCarlo=.TRUE., so mread applies the var file), all to
18 Gyr from the same var file."""
import gzip, io, os, sys, tarfile
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
# readers copied from cmp_tracks.py (which runs on import)
XCOLS = ['star_age','star_mass','star_mdot','he_core_mass','c_core_mass',
         'log_L','log_LH','log_LHe','log_Teff','log_R','log_g',
         'surface_h1','surface_he3','surface_he4','surface_c12','surface_o16',
         'log_center_T','log_center_Rho','center_gamma',
         'center_h1','center_he4','center_c12']

def read_dsepX(path_or_bytes):
    """dsepX writes MESA history ONE VALUE PER LINE."""
    if isinstance(path_or_bytes, bytes):
        txt = path_or_bytes.decode('utf-8', 'replace')
    else:
        op = gzip.open if str(path_or_bytes).endswith('.gz') else open
        with op(path_or_bytes, 'rt', errors='replace') as fh: txt = fh.read()
    lines = txt.split('\n')
    start = next((i+1 for i,l in enumerate(lines) if 'star_age' in l and 'log_Teff' in l), None)
    if start is None: return None
    vals = []
    for l in lines[start:]:
        s = l.strip()
        if not s: continue
        try: vals.append(float(s))
        except ValueError: continue
    n = len(vals)//len(XCOLS)
    if n == 0: return None
    a = np.array(vals[:n*len(XCOLS)]).reshape(n, len(XCOLS))
    return {c: a[:,j] for j,c in enumerate(XCOLS)}

def read_dsep4(path):
    """dsep4 writes native MESA: row per model."""
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt', errors='replace') as fh: lines = fh.read().split('\n')
    i = next((k for k,l in enumerate(lines) if 'star_age' in l and 'log_Teff' in l), None)
    if i is None: return None
    h = lines[i].split(); rows=[]
    for l in lines[i+1:]:
        p = l.split()
        if len(p) == len(h):
            try: rows.append([float(x) for x in p])
            except ValueError: pass
    if not rows: return None
    a = np.array(rows)
    return {c: a[:,j] for j,c in enumerate(h)}


RUN = "20000"
S = "/dartfs-hpc/scratch/f007058"
XT = f"{S}/Tracks_noco/mcfeh031.{RUN}.tar"
D4 = f"{S}/Tracks_dsep4_parts_18Gyr/{RUN}"
D4MC = f"{S}/Tracks_dsep4_mc_18Gyr/{RUN}"
OUT = f"{S}/dsep4_batch/compare"
MASSES = ["0700", "0800", "0850", "1000", "1200", "1500", "2000", "3000"]

xmem = {}
with tarfile.open(XT) as tf:
    for m in tf.getmembers():
        b = os.path.basename(m.name)
        if b.startswith("m") and b.endswith(".gz"):
            xmem[b] = gzip.decompress(tf.extractfile(m).read())


def load(M):
    key = f"m{M}.feh031.{RUN}.gz"
    tx = read_dsepX(xmem[key]) if key in xmem else None
    t4 = read_dsep4(f"{D4}/{key}") if os.path.exists(f"{D4}/{key}") else None
    tm = read_dsep4(f"{D4MC}/{key}") if os.path.exists(f"{D4MC}/{key}") else None
    return tx, t4, tm


def h_exhaust(t):
    i = np.argmax(t["center_h1"] < 1e-3)
    return t["star_age"][i] / 1e9 if t["center_h1"][i] < 1e-3 else np.nan


def rgb_teff(t, lo=1.0, hi=2.4):
    """log Teff on the giant branch at log L = lo..hi (after H exhaustion, ascending part)"""
    i0 = np.argmax(t["center_h1"] < 1e-3)
    L, T = t["log_L"][i0:], t["log_Teff"][i0:]
    k = np.flatnonzero((L > lo) & (L < hi))
    if len(k) < 5:
        return None
    L, T = L[k], T[k]
    o = np.argsort(L)
    grid = np.linspace(lo + 0.05, hi - 0.05, 20)
    return np.interp(grid, L[o], T[o])


fig, axes = plt.subplots(2, 4, figsize=(21, 10))
rows = []
for ax, M in zip(axes.ravel(), MASSES):
    tx, t4, tm = load(M)
    for t, lab, kw in [(tx, "dsepX_fixed", dict(color="#1f77b4", lw=1.6)),
                       (t4, "dsep4, MC off", dict(color="#d62728", lw=1.0, ls="--", alpha=0.7)),
                       (tm, "dsep4, MC on", dict(color="#2ca02c", lw=1.4))]:
        if t is not None:
            ax.plot(t["log_Teff"], t["log_L"], label=f"{lab} ({len(t['log_L'])} mod)", **kw)
    ax.set_title(f"{int(M)/1000:.2f} $M_\\odot$"); ax.invert_xaxis(); ax.grid(alpha=.25)
    ax.set_xlabel("log $T_{eff}$"); ax.set_ylabel("log $L/L_\\odot$"); ax.legend(fontsize=8)
    r = {"mass": int(M) / 1000}
    for t, k in [(tx, "X"), (t4, "4off"), (tm, "4on")]:
        if t is None:
            continue
        r[f"maxL_{k}"] = t["log_L"].max(); r[f"end_{k}"] = t["star_age"][-1] / 1e9
        r[f"Hex_{k}"] = h_exhaust(t); r[f"n_{k}"] = len(t["log_L"])
    gx = rgb_teff(tx) if tx is not None else None
    for t, k in [(t4, "4off"), (tm, "4on")]:
        g = rgb_teff(t) if t is not None else None
        if gx is not None and g is not None:
            r[f"dTeffRGB_{k}"] = np.median(g - gx)
    rows.append(r)
fig.suptitle(f"run {RUN}: dsepX_fixed vs dsep4 with Monte Carlo mode off and on (same var file, endage 18 Gyr)", fontsize=14)
fig.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(f"{OUT}/tracks_{RUN}_dsepX_vs_dsep4_mc.png", dpi=115)
print(f"saved {OUT}/tracks_{RUN}_dsepX_vs_dsep4_mc.png\n")
print(f"{'mass':>5} | {'max log L':^20} | {'end age (Gyr)':^23} | {'H exhaustion (Gyr)':^23} | {'RGB dlogTeff vs X':^17}")
print(f"{'':>5} | {'X':>6}{'4off':>7}{'4on':>7} | {'X':>7}{'4off':>8}{'4on':>8} | {'X':>7}{'4off':>8}{'4on':>8} | {'4off':>8}{'4on':>9}")
for r in rows:
    g = lambda k, f="{:7.3f}": (f.format(r[k]) if k in r and np.isfinite(r[k]) else "    -  ")
    print(f"{r['mass']:5.2f} | {g('maxL_X','{:6.3f}')}{g('maxL_4off')}{g('maxL_4on')} | {g('end_X')}{g('end_4off','{:8.3f}')}{g('end_4on','{:8.3f}')} | "
          f"{g('Hex_X')}{g('Hex_4off','{:8.3f}')}{g('Hex_4on','{:8.3f}')} | {g('dTeffRGB_4off','{:+8.4f}')}{g('dTeffRGB_4on','{:+9.4f}')}")
