#!/usr/bin/env python3
"""NGC 6791 isochrones: OG vs New vs Emulated."""
import os, glob
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

C = ["EEP","isochrone_age_yr","initial_mass","star_mass","star_mdot","he_core_mass",
     "c_core_mass","log_L","log_LH","log_LHe","log_Teff","log_R","log_g","surface_h1",
     "surface_he3","surface_he4","surface_c12","surface_o16","log_center_T",
     "log_center_Rho","center_gamma","center_h1","center_he4","center_c12"]
SRC = {
 'OG':       ('/dartfs-hpc/rc/lab/C/ChaboyerB/Isochrones_Only_July_10/{r}/isos/isochrones_example.txt', '#1f77b4'),
 'New':      ('/dartfs-hpc/scratch/f007058/Isochrones_noco/{r}/isos/isochrones_example.txt',            '#2ca02c'),
 'Emulated': ('/dartfs-hpc/scratch/f007058/Isochrones_emul4k_tip/{r}/isos/isochrones_example.txt',      '#d62728'),
}
OUT = '/dartfs-hpc/scratch/f007058/iso_plots'
WLO, WHI = 0.073, 0.761          # the sCMD fit window in log L

def load(p):
    d = pd.read_csv(p, sep=r'\s+', comment='#', skiprows=10); d.columns = C; return d

def slice_age(d, age):
    a = np.unique(d.isochrone_age_yr)
    return d[d.isochrone_age_yr == a[np.argmin(abs(a - age*1e9))]].sort_values('initial_mass')

def runs_with_all(n):
    out = []
    for r in range(20000, 21000):
        if all(os.path.exists(t.format(r=r)) for t, _ in SRC.values()):
            out.append(str(r))
        if len(out) >= n: break
    return out

# ---------- Figure 1: grid of runs at 8.5 Gyr, full HRD ----------
rs = runs_with_all(12)
print(f"runs with all three: {len(rs)} -> {rs[:4]}...")
for tag, (ylim, xlim, fn) in {
        'full':  ((None, None), (None, None), 'fig1_3way_8.5Gyr_full.png'),
        'window':((-0.1, 1.0), (3.82, 3.66),  'fig2_3way_8.5Gyr_window.png')}.items():
    fig, axes = plt.subplots(3, 4, figsize=(21, 14))
    for ax, r in zip(axes.ravel(), rs):
        for lab, (tpl, col) in SRC.items():
            try: s = slice_age(load(tpl.format(r=r)), 8.5)
            except Exception: continue
            ax.plot(s.log_Teff, s.log_L, '-', color=col, lw=1.4, label=lab)
        ax.set_title(f'run {r}', fontsize=10); ax.grid(alpha=.25)
        if tag == 'full': ax.invert_xaxis()
        else:
            ax.set_xlim(*xlim); ax.set_ylim(*ylim)
            ax.axhspan(WLO, WHI, color='gold', alpha=.12, zorder=0)
        ax.set_xlabel('log $T_{eff}$', fontsize=8); ax.set_ylabel('log $L/L_\\odot$', fontsize=8)
        ax.tick_params(labelsize=8)
    axes.ravel()[0].legend(fontsize=9)
    ttl = 'full HRD' if tag == 'full' else 'sCMD fit window (shaded)'
    fig.suptitle(f'NGC 6791 isochrones at 8.5 Gyr — OG vs New vs Emulated — {ttl}', fontsize=15)
    fig.tight_layout(rect=[0, 0, 1, 0.965]); fig.savefig(f'{OUT}/{fn}', dpi=105); plt.close(fig)
    print('saved', fn)

# ---------- Figure 3: one run across ages ----------
AGES = [4.0, 6.0, 8.0, 9.0, 10.0, 12.0]
fig, axes = plt.subplots(2, 3, figsize=(19, 11))
for ax, age in zip(axes.ravel(), AGES):
    for lab, (tpl, col) in SRC.items():
        try: s = slice_age(load(tpl.format(r=rs[0])), age)
        except Exception: continue
        ax.plot(s.log_Teff, s.log_L, '-', color=col, lw=1.5,
                label=f"{lab} (m≤{s.initial_mass.max():.2f})")
    ax.set_title(f'{age:.1f} Gyr'); ax.invert_xaxis(); ax.grid(alpha=.25)
    ax.set_xlabel('log $T_{eff}$'); ax.set_ylabel('log $L/L_\\odot$'); ax.legend(fontsize=8)
fig.suptitle(f'NGC 6791 run {rs[0]} — isochrone vs age — OG vs New vs Emulated', fontsize=15)
fig.tight_layout(rect=[0, 0, 1, 0.96]); fig.savefig(f'{OUT}/fig3_3way_ages.png', dpi=105); plt.close(fig)
print('saved fig3_3way_ages.png')

# ---------- Figure 4: quantitative offsets in the fit window ----------
many = runs_with_all(60)
rows = []
grid = np.linspace(WLO, WHI, 80)
for r in many:
    try:
        cur = {}
        for lab, (tpl, _) in SRC.items():
            s = slice_age(load(tpl.format(r=r)), 8.5).drop_duplicates('log_L').sort_values('log_L')
            cur[lab] = np.interp(grid, s.log_L, s.log_Teff, left=np.nan, right=np.nan)
        if any(np.isnan(v).all() for v in cur.values()): continue
        rows.append({'run': r,
                     'New-OG':   np.nanmean(cur['New'] - cur['OG']),
                     'Emul-OG':  np.nanmean(cur['Emulated'] - cur['OG']),
                     'Emul-New': np.nanmean(cur['Emulated'] - cur['New'])})
    except Exception: pass
df = pd.DataFrame(rows)
fig, ax = plt.subplots(1, 2, figsize=(14, 5.5))
for k, c in [('New-OG', '#2ca02c'), ('Emul-OG', '#d62728'), ('Emul-New', '#9467bd')]:
    ax[0].hist(df[k], bins=24, histtype='step', lw=2, color=c,
               label=f"{k}: {df[k].mean():+.4f} ± {df[k].std():.4f} dex")
ax[0].axvline(0, color='k', lw=.8, ls=':')
ax[0].set_xlabel('mean Δ log $T_{eff}$ in the fit window'); ax[0].set_ylabel('runs')
ax[0].set_title(f'Teff offset at 8.5 Gyr ({len(df)} runs)'); ax[0].legend(fontsize=8); ax[0].grid(alpha=.25)
for k, c in [('New-OG', '#2ca02c'), ('Emul-OG', '#d62728'), ('Emul-New', '#9467bd')]:
    ax[1].hist(5500*(10**df[k]-1), bins=24, histtype='step', lw=2, color=c,
               label=f"{k}: {5500*(10**df[k].mean()-1):+.0f} K")
ax[1].axvline(0, color='k', lw=.8, ls=':')
ax[1].set_xlabel('ΔT$_{eff}$ at 5500 K  [K]'); ax[1].set_ylabel('runs')
ax[1].set_title('same, in kelvin'); ax[1].legend(fontsize=8); ax[1].grid(alpha=.25)
fig.tight_layout(); fig.savefig(f'{OUT}/fig4_offsets.png', dpi=115); plt.close(fig)
print('saved fig4_offsets.png')
print(df[['New-OG','Emul-OG','Emul-New']].describe().loc[['mean','std']].round(5).to_string())
df.to_csv(f'{OUT}/offsets_8.5Gyr.csv', index=False)
