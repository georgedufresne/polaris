#!/usr/bin/env python3
"""DSEE stage-2 training tracks with EEP anchors, 10 tracks per figure.

Each tarball is one random (physics, mass) draw, so a figure is 10 independent
stars rather than a mass sequence -- the label gives each one's mass.
"""
import gzip, os, sys, tarfile
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# reuse the EXACT detector from plot_eeps.py without running its plotting body
_src = open('/dartfs-hpc/scratch/f007058/dsep4_batch/compare/plot_eeps.py').read()
_ns = {}; exec(_src[:_src.index('files = sorted(')], _ns)
find_eeps = _ns['find_eeps']; EEPS = _ns['EEPS']

D   = '/dartfs-hpc/scratch/f007058/stage_2_tracks'
OUT = '/dartfs-hpc/scratch/f007058/dsee_eep_plots'
NPLOT   = int(sys.argv[1]) if len(sys.argv) > 1 else 30     # number of figures
PERPLOT = int(sys.argv[2]) if len(sys.argv) > 2 else 10     # tracks per figure
os.makedirs(OUT, exist_ok=True)

def read_mesa(buf):
    lines = buf.decode('utf-8', 'replace').split('\n')
    i = next((k for k, l in enumerate(lines) if 'star_age' in l and 'log_Teff' in l), None)
    if i is None: return None
    h = lines[i].split(); rows = []
    for l in lines[i+1:]:
        q = l.split()
        if len(q) == len(h):
            try: rows.append([float(x) for x in q])
            except ValueError: pass
    if len(rows) < 20: return None
    a = np.array(rows)
    return {c: a[:, j] for j, c in enumerate(h)}

def load(tarname):
    with tarfile.open(f'{D}/{tarname}') as tf:
        m = next((x for x in tf.getmembers()
                  if os.path.basename(x.name).startswith('m')
                  and not os.path.basename(x.name).startswith('iso_')), None)
        if m is None: return None, None
        t = read_mesa(gzip.decompress(tf.extractfile(m).read()))
    if t is None: return None, None
    return float(os.path.basename(m.name)[1:5]) / 1000.0, t

# os.listdir() materialises the whole directory and stalls on this NFS share;
# scandir yields lazily so we can stop as soon as we have enough candidates.
need = NPLOT * PERPLOT
CAP = need * 4
tars = []
with os.scandir(D) as it:
    for ent in it:
        if ent.name.endswith('.tar'):
            tars.append(ent.name)
            if len(tars) >= CAP: break
rng = np.random.default_rng(1)
rng.shuffle(tars)
pick = tars
print(f"scanned {len(tars)} tarballs (capped); need {need} readable tracks", flush=True)

COL = plt.cm.tab10.colors
TALLY = {k: 0 for k, _, _ in EEPS}; STATUS = {}; WHICH = {}; HEC = []
buf, fig_i = [], 0
for f in pick:
    if fig_i >= NPLOT: break
    try: mass, t = load(f)
    except Exception: continue
    if t is None: continue
    buf.append((f, mass, t))
    if len(buf) < PERPLOT: continue

    fig, ax = plt.subplots(figsize=(10.5, 8.5))
    for j, (fn, mm, tt) in enumerate(buf):
        e = find_eeps(tt)
        st = e.get('_status', 'truncated'); STATUS[st] = STATUS.get(st, 0) + 1
        if '_which' in e: WHICH[e['_which']] = WHICH.get(e['_which'], 0) + 1
        if '_hecore' in e: HEC.append(e['_hecore'])
        z = e.get('ZAMS', 0)
        c = COL[j % len(COL)]
        ax.plot(tt['log_Teff'][z:], tt['log_L'][z:], '-', color=c, lw=1.1, alpha=0.9,
                label=f"{fn.replace('DSEE.','').replace('.tar','')}  {mm:.3f} M$_\\odot$  [{st}]")
        for name, mk, mc in EEPS:
            if name.startswith('_') or name not in e: continue
            k = e[name]
            if k < z: continue
            ax.scatter(tt['log_Teff'][k], tt['log_L'][k], marker=mk,
                       s=135 if mk == '*' else 46, facecolor=mc,
                       edgecolor='white', lw=0.7, zorder=7)
            TALLY[name] += 1
    ax.invert_xaxis(); ax.grid(alpha=.25)
    ax.set_xlabel('log $T_{eff}$'); ax.set_ylabel('log $L/L_\\odot$')
    ax.set_title(f'DSEE stage-2 tracks, set {fig_i+1} of {NPLOT} — EEP anchors')
    l1 = ax.legend(fontsize=7, loc='lower left', title='track  (terminus status)')
    ax.add_artist(l1)
    ax.legend(handles=[Line2D([], [], marker=mk, color='w', markerfacecolor=mc,
                              markeredgecolor='k', markersize=11 if mk == '*' else 8, label=name)
                       for name, mk, mc in EEPS], fontsize=8, loc='upper right', title='EEP')
    fig.tight_layout()
    p = f'{OUT}/dsee_eeps_set{fig_i+1:02d}.png'
    fig.savefig(p, dpi=115); plt.close(fig)
    print(f'  saved {os.path.basename(p)}', flush=True)
    buf = []; fig_i += 1

print(f"\nfigures: {fig_i}   tracks: {fig_i*PERPLOT}")
print('  anchors: ' + '  '.join(f'{k}:{v}' for k, v in TALLY.items()))
print('  fired by: ' + '  '.join(f'{k[5:]}={v}' for k, v in sorted(WHICH.items())))
print('  status:   ' + '  '.join(f'{k}={v}' for k, v in sorted(STATUS.items())) +
      (f'   he_core at terminus: median {np.median(HEC):.3f}' if HEC else ''))
