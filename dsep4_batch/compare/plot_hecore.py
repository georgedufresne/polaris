#!/usr/bin/env python3
"""All 39 NGC 6791 tracks (0.20-3.00 Msun) on one HR diagram, with markers where
the helium core first reaches 0.1, 0.2 and 0.3 Msun."""
import gzip, os
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

RUN = '20000'
D = f'/dartfs-hpc/scratch/f007058/Tracks_dsep4_parts/{RUN}'
OUT = '/dartfs-hpc/scratch/f007058/dsep4_batch/compare'
LEVELS = [(0.1, 'o', 'He core 0.1 M$_\\odot$'),
          (0.2, 's', 'He core 0.2 M$_\\odot$'),
          (0.3, '^', 'He core 0.3 M$_\\odot$')]
MCOL = {0.1: '#111111', 0.2: '#d62728', 0.3: '#1f77b4'}

def read(p):
    with gzip.open(p, 'rt', errors='replace') as fh:
        lines = fh.read().split('\n')
    i = next((k for k, l in enumerate(lines) if 'star_age' in l and 'log_Teff' in l), None)
    if i is None: return None
    h = lines[i].split(); rows = []
    for l in lines[i+1:]:
        p2 = l.split()
        if len(p2) == len(h):
            try: rows.append([float(x) for x in p2])
            except ValueError: pass
    if len(rows) < 10: return None
    a = np.array(rows)
    return {c: a[:, j] for j, c in enumerate(h)}

files = sorted(f for f in os.listdir(D) if f.startswith('m') and f.endswith('.gz'))
cmap = plt.cm.viridis
masses = [int(f[1:5])/1000 for f in files]
norm = plt.Normalize(min(masses), max(masses))

fig, ax = plt.subplots(figsize=(11, 9))
hits = {0.1: 0, 0.2: 0, 0.3: 0}
for f in files:
    m = int(f[1:5]) / 1000
    t = read(f'{D}/{f}')
    if t is None: continue
    ax.plot(t['log_Teff'], t['log_L'], '-', color=cmap(norm(m)), lw=0.9, alpha=0.85)
    hc = t.get('he_core_mass')
    if hc is None: continue
    for lev, mk, _ in LEVELS:
        idx = np.nonzero(hc >= lev)[0]          # first model reaching this core mass
        if len(idx):
            k = idx[0]
            ax.scatter(t['log_Teff'][k], t['log_L'][k], marker=mk, s=58,
                       facecolor=MCOL[lev], edgecolor='white', lw=0.8, zorder=6)
            hits[lev] += 1

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
cb = fig.colorbar(sm, ax=ax, pad=0.015); cb.set_label('initial mass [M$_\\odot$]')
handles = [Line2D([], [], marker=mk, color='w', markerfacecolor=MCOL[lev],
                  markeredgecolor='k', markersize=9, label=f'{lab}  (n={hits[lev]})')
           for lev, mk, lab in LEVELS]
ax.legend(handles=handles, loc='lower left', fontsize=10)
ax.invert_xaxis(); ax.grid(alpha=.25)
ax.set_xlabel('log $T_{eff}$'); ax.set_ylabel('log $L/L_\\odot$')
ax.set_title(f'NGC 6791 run {RUN} (dsep4): 39 tracks, 0.20-3.00 M$_\\odot$\n'
             'markers = first model where the He core reaches 0.1 / 0.2 / 0.3 M$_\\odot$')
fig.tight_layout()
o = f'{OUT}/tracks_hecore_{RUN}.png'
fig.savefig(o, dpi=125); print('saved', o)
print(f"tracks plotted: {len(files)}   markers: 0.1 -> {hits[0.1]}, 0.2 -> {hits[0.2]}, 0.3 -> {hits[0.3]}")
