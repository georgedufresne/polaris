#!/usr/bin/env python3
"""EEP anchors on a sample of the DSEE training tracks (Catherine's grid -- what W was trained on).

Unlike the cluster grids, DSEE stores ONE mass per tarball: each MC set is a random
(physics, mass) draw, so this is a scatter of independent tracks across the whole
training domain rather than a mass sequence at fixed physics.
"""
import gzip, io, os, sys, tarfile
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# reuse the EXACT detector from plot_eeps.py without running its plotting body
_src = open('/dartfs-hpc/scratch/f007058/dsep4_batch/compare/plot_eeps.py').read()
_ns = {}
exec(_src[:_src.index('files = sorted(')], _ns)
find_eeps = _ns['find_eeps']; EEPS = _ns['EEPS']

D = os.environ.get('DSEE_DIR', '/dartfs-hpc/scratch/f007058/stage_2_tracks')
N = int(sys.argv[1]) if len(sys.argv) > 1 else 200
OUT = '/dartfs-hpc/scratch/f007058/dsep4_batch/compare'

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

tars = sorted(f for f in os.listdir(D) if f.endswith('.tar'))
print(f"DSEE tarballs available: {len(tars):,}; sampling {N}")
rng = np.random.default_rng(0)
pick = [tars[i] for i in rng.choice(len(tars), size=min(N, len(tars)), replace=False)]

tracks = []
for f in pick:
    try:
        with tarfile.open(f'{D}/{f}') as tf:
            m = next((x for x in tf.getmembers()
                      if os.path.basename(x.name).startswith('m')
                      and not os.path.basename(x.name).startswith('iso_')), None)
            if m is None: continue
            t = read_mesa(gzip.decompress(tf.extractfile(m).read()))
        if t is None: continue
        mass = float(os.path.basename(m.name)[1:5]) / 1000.0
        tracks.append((mass, t))
    except Exception:
        continue
print(f"readable tracks: {len(tracks)}")
if not tracks:
    sys.exit(0)

masses = [m for m, _ in tracks]
cmap = plt.cm.viridis; norm = plt.Normalize(min(masses), max(masses))
fig, ax = plt.subplots(figsize=(12, 9.5))
counts = {k: 0 for k, _, _ in EEPS}; STATUS = {}; HEC = []; WHICH = {}
for mass, t in tracks:
    e = find_eeps(t)
    st = e.get('_status', 'truncated'); STATUS[st] = STATUS.get(st, 0) + 1
    if '_hecore' in e: HEC.append(e['_hecore'])
    if '_which' in e: WHICH[e['_which']] = WHICH.get(e['_which'], 0) + 1
    z = e.get('ZAMS', 0)
    ax.plot(t['log_Teff'][z:], t['log_L'][z:], '-', color=cmap(norm(mass)), lw=0.5, alpha=0.5)
    for name, mk, col in EEPS:
        if name.startswith('_') or name not in e: continue
        k = e[name]
        if k < z: continue
        ax.scatter(t['log_Teff'][k], t['log_L'][k], marker=mk,
                   s=90 if mk == '*' else 24, facecolor=col,
                   edgecolor='white', lw=0.4, zorder=6)
        counts[name] += 1

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
cb = fig.colorbar(sm, ax=ax, pad=0.015); cb.set_label('initial mass [M$_\\odot$]')
handles = [Line2D([], [], marker=mk, color='w', markerfacecolor=col, markeredgecolor='k',
                  markersize=12 if mk == '*' else 8, label=f'{name}  (n={counts[name]})')
           for name, mk, col in EEPS]
ax.legend(handles=handles, loc='lower left', fontsize=10, title='EEP anchors')
ax.invert_xaxis(); ax.grid(alpha=.25)
ax.set_xlabel('log $T_{eff}$'); ax.set_ylabel('log $L/L_\\odot$')
ax.set_title(f'DSEE stage-2 training tracks\n'
             f'{len(tracks)} random (physics, mass) draws — EEP anchors')
fig.tight_layout()
o = f'{OUT}/dsee_eeps_stage2_{len(tracks)}.png'
fig.savefig(o, dpi=120); print('saved', o)
print('  ' + '  '.join(f'{k}:{v}' for k, v in counts.items()))
print('  fired by: ' + '  '.join(f'{k[5:]}={v}' for k, v in sorted(WHICH.items())))
print('  status:   ' + '  '.join(f'{k}={v}' for k, v in sorted(STATUS.items())) +
      (f'   he_core at terminus: median {np.median(HEC):.3f}' if HEC else ''))
print(f'  mass range sampled: {min(masses):.3f} - {max(masses):.3f}')
