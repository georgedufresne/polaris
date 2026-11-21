#!/usr/bin/env python3
"""dsepX_fixed vs dsep4 stellar tracks -- same MC parameters, same endage (18 Gyr).

Isolates the CODE difference: both grids use varfeh+031.<run> with identical
physics multipliers and endage=1.8D10, so any difference is dsepX vs dsep4.
"""
import gzip, os, sys, tarfile, tempfile
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

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

RUN = sys.argv[1] if len(sys.argv) > 1 else '20000'
MASSES = ['0700','0800','0850','1000','1200','1500','2000','3000']
D4 = f'/dartfs-hpc/scratch/f007058/Tracks_dsep4_parts_18Gyr/{RUN}'
XT = f'/dartfs-hpc/scratch/f007058/Tracks_noco/mcfeh031.{RUN}.tar'

# pull dsepX tracks out of the tar
xmem = {}
if os.path.exists(XT):
    with tarfile.open(XT) as tf:
        for m in tf.getmembers():
            b = os.path.basename(m.name)
            if b.startswith('m') and b.endswith('.gz'):
                xmem[b] = gzip.decompress(tf.extractfile(m).read())
print(f"run {RUN}: dsepX tracks in tar = {len(xmem)}, dsep4 files = "
      f"{len(os.listdir(D4)) if os.path.isdir(D4) else 0}")

fig, axes = plt.subplots(2, 4, figsize=(21, 10), sharex=False)
stats = []
for ax, M in zip(axes.ravel(), MASSES):
    key = f'm{M}.feh031.{RUN}.gz'
    tx = read_dsepX(xmem[key]) if key in xmem else None
    p4 = f'{D4}/m{M}.feh031.{RUN}.gz'
    t4 = read_dsep4(p4) if os.path.exists(p4) else None
    if tx is not None:
        ax.plot(tx['log_Teff'], tx['log_L'], '-', lw=1.4, color='#1f77b4',
                label=f"dsepX ({len(tx['log_L'])} mod)")
    if t4 is not None:
        ax.plot(t4['log_Teff'], t4['log_L'], '-', lw=1.4, color='#d62728',
                label=f"dsep4 ({len(t4['log_L'])} mod)")
    ax.set_title(f"{int(M)/1000:.2f} $M_\\odot$", fontsize=11)
    ax.invert_xaxis(); ax.grid(alpha=.25)
    ax.set_xlabel('log $T_{eff}$'); ax.set_ylabel('log $L/L_\\odot$')
    ax.legend(fontsize=8, loc='best')
    if tx is not None and t4 is not None:
        stats.append((M, tx['log_L'].max(), t4['log_L'].max(),
                      tx['star_age'][-1]/1e9, t4['star_age'][-1]/1e9,
                      len(tx['log_L']), len(t4['log_L'])))
fig.suptitle(f'NGC 6791 tracks, run {RUN}: dsepX_fixed vs dsep4 '
             f'(identical MC parameters, endage 18 Gyr)', fontsize=14)
fig.tight_layout(rect=[0,0,1,0.96])
out = f'/dartfs-hpc/scratch/f007058/dsep4_batch/compare/tracks_{RUN}_dsepX_vs_dsep4.png'
fig.savefig(out, dpi=115); print(f"saved {out}")

# --- turnoff zoom: the region the sCMD fit actually uses (log L 0.07-0.76) ---
fig2, ax2 = plt.subplots(1, 2, figsize=(14, 6))
for M, col in zip(['0850','1000','1200'], ['#1f77b4','#2ca02c','#9467bd']):
    key=f'm{M}.feh031.{RUN}.gz'
    tx = read_dsepX(xmem[key]) if key in xmem else None
    p4 = f'{D4}/m{M}.feh031.{RUN}.gz'
    t4 = read_dsep4(p4) if os.path.exists(p4) else None
    if tx is not None: ax2[0].plot(tx['log_Teff'], tx['log_L'], '-',  color=col, lw=1.5, label=f'{int(M)/1000:.2f} M$_\\odot$')
    if t4 is not None: ax2[0].plot(t4['log_Teff'], t4['log_L'], '--', color=col, lw=1.5)
    if tx is not None: ax2[1].plot(tx['log_Teff'], tx['log_L'], '-',  color=col, lw=1.8, label=f'{int(M)/1000:.2f} dsepX')
    if t4 is not None: ax2[1].plot(t4['log_Teff'], t4['log_L'], '--', color=col, lw=1.8, label=f'{int(M)/1000:.2f} dsep4')
ax2[0].set_title('full tracks (solid = dsepX, dashed = dsep4)')
ax2[1].set_title('zoom: the sCMD fit window, log L 0.07-0.76')
ax2[1].set_ylim(0.0, 0.85); ax2[1].set_xlim(3.80, 3.68)
for a in ax2:
    a.invert_xaxis() if a is ax2[0] else None
    a.grid(alpha=.25); a.set_xlabel('log $T_{eff}$'); a.set_ylabel('log $L/L_\\odot$')
    a.legend(fontsize=8)
ax2[1].axhspan(0.073, 0.761, color='gold', alpha=.12, zorder=0)
fig2.suptitle(f'run {RUN}: where dsepX and dsep4 agree and disagree', fontsize=13)
fig2.tight_layout(rect=[0,0,1,0.94])
o2=f'/dartfs-hpc/scratch/f007058/dsep4_batch/compare/tracks_{RUN}_turnoff_zoom.png'
fig2.savefig(o2, dpi=115); print(f"saved {o2}")

print(f"\n{'mass':>6}{'dsepX maxL':>12}{'dsep4 maxL':>12}{'dL':>8}"
      f"{'dsepX age':>11}{'dsep4 age':>11}{'dsepX n':>9}{'dsep4 n':>9}")
for M,lx,l4,ax_,a4,nx,n4 in stats:
    print(f"{int(M)/1000:>6.2f}{lx:>12.3f}{l4:>12.3f}{l4-lx:>+8.3f}"
          f"{ax_:>11.3f}{a4:>11.3f}{nx:>9d}{n4:>9d}")
