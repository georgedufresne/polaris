#!/usr/bin/env python3
"""All 39 NGC 6791 tracks on one HR diagram, with the six v1.2 primary EEPs marked.

  1 PreMS    log Tc >= 5.0
  2 ZAMS     max log g before Xc has dropped 0.001
  3 IAMS     Xc = 0.30
  4 TAMS     Xc = 1e-12
  5 RGBBase  the HR knee between TAMS and the tip
  6 RGBTip   later of max log L / min log Teff, HeIgn = log(L_He/L) >= -2
"""
import gzip, os, sys, tarfile
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

RUN = sys.argv[1] if len(sys.argv) > 1 else '10000'
GATE = -3.5          # log(L_He/L_H) threshold for tagging an RGB tip (was -3.0)
TAR = f'/dartfs-hpc/scratch/f007058/Tracks_tuc_noco/mcfeh077.{RUN}.tar'
OUT = '/dartfs-hpc/scratch/f007058/dsep4_batch/compare'

EEPS = [('ZAMS',       's', '#1f77b4'),
        ('TAMS',       '^', '#ff7f0e'),
        ('Xc0.95'.replace('Xc0.','Xc.'), '.', '#0d0887'),
        ('Xc0.9'.replace('Xc0.','Xc.'), '.', '#3a049a'),
        ('Xc0.8'.replace('Xc0.','Xc.'), '.', '#5c01a6'),
        ('Xc0.7'.replace('Xc0.','Xc.'), '.', '#7e03a8'),
        ('Xc0.6'.replace('Xc0.','Xc.'), '.', '#9c179e'),
        ('Xc0.5'.replace('Xc0.','Xc.'), '.', '#b52f8c'),
        ('Xc0.4'.replace('Xc0.','Xc.'), '.', '#cc4778'),
        ('Xc0.3'.replace('Xc0.','Xc.'), '.', '#de5f65'),
        ('Xc0.2'.replace('Xc0.','Xc.'), '.', '#ed7953'),
        ('Xc0.1'.replace('Xc0.','Xc.'), '.', '#f89540'),
        ('Xc0.01'.replace('Xc0.','Xc.'), '.', '#fdb42f'),
        ('Xc0.001'.replace('Xc0.','Xc.'), '.', '#fbd524'),
        ('Xc1e-06'.replace('Xc0.','Xc.'), '.', '#f0f921'),
        ('HeL_m10',    'o', '#0000ff'),
        ('HeL_m8',    'o', '#0040df'),
        ('HeL_m6',    'o', '#0080bf'),
        ('HeL_m4',    'o', '#00c09f'),
        ('HeL_m2',    'o', '#00ff80'),
        ('RGBBase_s',  '>', '#8c564b')]   # slope threshold, kept only if HeL_m10 is reached

XC = ['star_age','star_mass','star_mdot','he_core_mass','c_core_mass','log_L','log_LH',
      'log_LHe','log_Teff','log_R','log_g','surface_h1','surface_he3','surface_he4',
      'surface_c12','surface_o16','log_center_T','log_center_Rho','center_gamma',
      'center_h1','center_he4','center_c12']

def read(buf):
    """dsepX writes MESA history ONE VALUE PER LINE."""
    lines = buf.decode('utf-8', 'replace').split('\n')
    i = next((k + 1 for k, l in enumerate(lines) if 'star_age' in l and 'log_Teff' in l), None)
    if i is None: return None
    v = []
    for l in lines[i:]:
        q = l.strip()
        if not q: continue
        try: v.append(float(q))
        except ValueError: continue
    n = len(v) // len(XC)
    if n < 20: return None
    a = np.array(v[:n * len(XC)]).reshape(n, len(XC))
    return {c: a[:, j] for j, c in enumerate(XC)}

def find_eeps(t):
    """Return {name: index} for whichever anchors exist on this track."""
    out = {}
    xc = t.get('center_h1'); lg = t.get('log_g'); lct = t.get('log_center_T')
    lL = t['log_L']; lT = t['log_Teff']
    lHe = t.get('log_LHe'); lH = t.get('log_LH')
    n = len(lL)

    if lct is not None:                                   # 1 PreMS
        k = np.nonzero(lct >= 5.0)[0]
        if len(k): out['PreMS'] = int(k[0])

    if xc is not None and lg is not None:                 # 2 ZAMS
        drop = np.nonzero(xc <= xc[0] - 0.001)[0]
        hi = int(drop[0]) if len(drop) else n
        if hi > 1: out['ZAMS'] = int(np.argmax(lg[:hi]))

    # --- main-sequence ladder: Xc / Xc(ZAMS) ---
    zi = out.get('ZAMS')
    if xc is not None and zi is not None and xc[zi] > 0:
        xz = xc[zi]
        for f in (0.95,0.90,0.80,0.70,0.60,0.50,0.40,0.30,0.20,0.10,0.01,1e-3,1e-6):
            k = np.nonzero(xc[zi:] <= f * xz)[0]
            if len(k): out[('Xc%g' % f).replace('Xc0.','Xc.')] = zi + int(k[0])

    if xc is not None:
        k = np.nonzero(xc <= 1e-12)[0]                    # 4 TAMS
        if len(k): out['TAMS'] = int(k[0])

    tams = out.get('TAMS')

    # --- log(L_He/L) >= -2 : helium supplies 1% of the total luminosity ---
    if tams is not None and lHe is not None and tams < n - 2:
        for lv, nm in ((-2.0,'HeL_m2'),(-4.0,'HeL_m4'),(-6.0,'HeL_m6'),(-8.0,'HeL_m8'),(-10.0,'HeL_m10')):
            k = np.nonzero((lHe[tams:] - lL[tams:]) >= lv)[0]
            if len(k): out[nm] = tams + int(k[0])

    # --- TERMINUS: four gated candidates, take the LATEST valid one. ---
    # L_max / Teff_min are only trusted if interior (>=5 rows from the end), so a
    # track that simply stops does not get a spurious extremum at its last row.
    # M_He and Tc are core criteria that still fire on the degenerate branch.
    if tams is not None and tams < n - 2:
        cand = {}
        if n - 1 - (tams + int(np.argmax(lL[tams:]))) >= 5:
            cand['Term_Lmax'] = tams + int(np.argmax(lL[tams:]))
        if n - 1 - (tams + int(np.argmin(lT[tams:]))) >= 5:
            cand['Term_Tmin'] = tams + int(np.argmin(lT[tams:]))
        hcm = t.get('he_core_mass')
        if hcm is not None:
            k = np.nonzero(hcm[tams:] >= 0.45)[0]
            if len(k): cand['Term_MHe'] = tams + int(k[0])
        if lct is not None:
            k = np.nonzero(lct[tams:] >= 7.80)[0]
            if len(k): cand['Term_Tc'] = tams + int(k[0])
        if cand:
            which = max(cand, key=lambda q: cand[q])        # kept for bookkeeping only
            out['_which'] = which
            out['_status'] = ('complete' if which in ('Term_Lmax', 'Term_Tmin')
                              else 'near_tip')
            out['_hecore'] = float(t['he_core_mass'][cand[which]]) if 'he_core_mass' in t else float('nan')
        else:
            out['_status'] = 'truncated'

    # --- RGBBase: first row after TAMS where L_H/L_total recovers above 0.95,
    # taken AFTER its post-TAMS minimum (H burning dips in the shell, then re-establishes).
    if tams is not None and lH is not None and n - tams > 5:
        r = 10.0 ** (lH[tams:] - lL[tams:])
        lo = int(np.argmin(r))                            # the post-TAMS minimum
        k = np.nonzero(r[lo:] > 0.95)[0]
        if len(k): out['RGBBase_LH'] = tams + lo + int(k[0])

    tip = out.get('RGBTip')
    # --- absolute central density ladder (log rho_c) ---
    lcr = t.get('log_center_Rho')
    if lcr is not None:
        for lv in (4.5, 5.5):
            k = np.nonzero(lcr >= lv)[0]
            if len(k): out[f'logRc{lv:.1f}'] = int(k[0])

    # --- Delta log g from TAMS: post-TAMS stars only expand, so log g falls
    # monotonically at every mass. Surface property -> no core-physics regime split,
    # and the range is large (~4 -> 0.5).
    if tams is not None and lg is not None and tams < n - 2:
        drop = lg[tams] - lg[tams:]
        for lv in (1.0, 2.0, 3.0):
            k = np.nonzero(drop >= lv)[0]
            if len(k): out[f'dlogg{lv:.1f}'] = tams + int(k[0])

    # --- RGBBase variant 2: maximum LOCAL curvature along TAMS -> tip ---
    # arc length in the distance metric (pyeep v1.2 uses Teff_scale = logL_scale = 2.0,
    # i.e. log L weighted like log Teff), theta = atan2(dlogL, -dlogTeff), take argmax|dtheta/ds|
    # over the FIRST HALF of the segment so the tip's own turn cannot win.
    if tams is not None and tip is not None and tip - tams > 12:
        TS, LS = 2.0, 2.0
        x = lT[tams:tip+1] * TS; y = lL[tams:tip+1] * LS
        dx = np.diff(x); dy = np.diff(y)
        ds = np.hypot(dx, dy)
        good = ds > 0
        if good.sum() > 10:
            th = np.unwrap(np.arctan2(dy[good], -dx[good]))
            ss = np.cumsum(ds[good])
            w = 5                                        # smooth; raw dtheta/ds is very noisy
            if len(th) > 3 * w:
                k = np.abs(np.gradient(th, ss))
                ker = np.ones(w) / w
                ks = np.convolve(k, ker, mode='same')
                half = len(ks) // 2
                idx_good = np.nonzero(good)[0]
                out['RGBBase_k'] = tams + int(idx_good[int(np.argmax(ks[:half]))])

    # --- RGBBase variant 3: first row where the track turns steeply vertical ---
    # dlogL/dlogTeff <= -8, smoothed over a few rows. Local, one tunable number, needs no tip.
    if tams is not None and n - tams > 12:
        hi = n - 1                                   # tip-independent by design: search all post-TAMS
        xs = lT[tams:hi+1]; ys = lL[tams:hi+1]
        if len(xs) > 12:
            w = 5
            ker = np.ones(w) / w
            xsm = np.convolve(xs, ker, mode='same'); ysm = np.convolve(ys, ker, mode='same')
            d = np.gradient(xsm)
            slope = np.where(np.abs(d) > 1e-9, np.gradient(ysm) / np.where(np.abs(d) > 1e-9, d, 1), 0.0)
            # PERSISTENCE: require the slope condition to hold for PERSIST consecutive
            # rows. The subgiant turn is transient; the turn onto the Hayashi line is not.
            PERSIST = 10
            cond = (slope[w:-w] <= -8.0).astype(int)
            if len(cond) > PERSIST:
                run = np.convolve(cond, np.ones(PERSIST, int), mode='valid')
                v = np.nonzero(run == PERSIST)[0]
                if len(v): out['RGBBase_s'] = tams + w + int(v[0])

    # RGBBase_s only counts as an EEP if the track actually gets onto the giant
    # branch, i.e. reaches log(L_He/L) >= -10. Without that gate it fires on the
    # subgiant turn of tracks that never become giants.
    if 'RGBBase_s' in out and 'HeL_m10' not in out:
        del out['RGBBase_s']

    if tams is not None and tip is not None and tip - tams > 5:   # 5 RGBBase: the knee
        x = lT[tams:tip+1]; y = lL[tams:tip+1]
        x0, y0, x1, y1 = x[0], y[0], x[-1], y[-1]
        dx, dy = x1 - x0, y1 - y0
        L = np.hypot(dx, dy)
        if L > 0:                                         # max perpendicular offset from the TAMS-tip chord
            d = np.abs(dy * (x - x0) - dx * (y - y0)) / L
            out['RGBBase'] = tams + int(np.argmax(d))
    return out

with tarfile.open(TAR) as tf:
    MEM = {os.path.basename(m.name): gzip.decompress(tf.extractfile(m).read())
           for m in tf.getmembers()
           if os.path.basename(m.name).startswith('m') and os.path.basename(m.name).endswith('.gz')}
files = sorted(MEM)
masses = [int(f[1:5]) / 1000 for f in files]
if not masses:
    print(f'  {RUN}: no readable tracks'); sys.exit(0)
cmap = plt.cm.viridis; norm = plt.Normalize(min(masses), max(masses))

fig, ax = plt.subplots(figsize=(12, 9.5))
counts = {k: 0 for k, _, _ in EEPS}
STATUS = {}; HEC = []; WHICH = {}
for f in files:
    m = int(f[1:5]) / 1000
    t = read(MEM[f])
    if t is None: continue
    e = find_eeps(t)
    STATUS[e.get('_status', 'truncated')] = STATUS.get(e.get('_status', 'truncated'), 0) + 1
    if '_hecore' in e: HEC.append(e['_hecore'])
    if '_which' in e: WHICH[e['_which']] = WHICH.get(e['_which'], 0) + 1
    z = e.get('ZAMS', 0)                      # start the curve at the ZAMS; drop PreMS->ZAMS
    ax.plot(t['log_Teff'][z:], t['log_L'][z:], '-', color=cmap(norm(m)), lw=0.9, alpha=0.85)
    for name, mk, col in EEPS:
        if name.startswith('_'): continue
        if name in e:
            k = e[name]
            if k < z: continue
            ax.scatter(t['log_Teff'][k], t['log_L'][k], marker=mk,
                       s=150 if mk == '*' else (14 if mk == '.' else 46), facecolor=col,
                       edgecolor='white', lw=0.7, zorder=7)
            counts[name] += 1

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm); sm.set_array([])
cb = fig.colorbar(sm, ax=ax, pad=0.015); cb.set_label('initial mass [M$_\\odot$]')
handles = [Line2D([], [], marker=mk, color='w', markerfacecolor=col, markeredgecolor='k',
                  markersize=13 if mk == '*' else 9, label=f'{name}  (n={counts[name]})')
           for i, (name, mk, col) in enumerate(EEPS)]
ax.legend(handles=handles, loc='lower left', fontsize=10, title='primary EEPs (pyeep v1.2)')
ax.invert_xaxis(); ax.grid(alpha=.25)
ax.set_xlabel('log $T_{eff}$'); ax.set_ylabel('log $L/L_\\odot$')
ax.set_title(f'47 Tuc run {RUN} (dsepX_fixed): 39 tracks, 0.20-3.00 M$_\\odot$ — from ZAMS onward, primary EEP anchors')
fig.tight_layout()
o = f'{OUT}/tuc_eeps_fromZAMS_{RUN}.png'
fig.savefig(o, dpi=125); print('saved', o)
print('  ' + '  '.join(f'{k}:{v}' for k, v in counts.items()) + f'   (of {len(files)} tracks)')
print('  fired by: ' + '  '.join(f'{k[5:]}={v}' for k, v in sorted(WHICH.items())))
print('  status: ' + '  '.join(f'{k}={v}' for k, v in sorted(STATUS.items())) +
      (f'   he_core at terminus: median {np.median(HEC):.3f}' if HEC else ''))
