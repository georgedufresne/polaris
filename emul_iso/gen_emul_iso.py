#!/usr/bin/env python3
"""Generate CENTRAL (sampling-free) emulated NGC 6791 isochrones from CONF1DENCE.

Central point: the conditional flow's base distribution is a standard Gaussian,
so its mode/mean is z = 0.  Pushing z = 0 through the inverse transform gives one
deterministic point per (physics, age, mass) label -- no RNG, no medoid, no
temperature.  This is `_sample_with_temperature(dist, 1, 0.0)` in
scripts/emulator_compare/iso_compare_m15.py, reduced to a single pass.

Output is MIST-format, byte-compatible with what the sCMD pipeline reads
(pd.read_csv(..., sep=r'\s+', comment='#', skiprows=10), 24 columns).
"""
import argparse, json, os, sys
import numpy as np
import pandas as pd
import torch, zuko

import CONF1DENCE
from CONF1DENCE.utils.helper import data_transform as dt

REPO = CONF1DENCE.__abspath__
CKPT = os.path.join(REPO, 'data/model/DSEE_run2_lr1e4/best.pt')
FMT  = os.path.join(REPO, 'data/format/DSEE_iso_format.json')

LOG_T_SUN, LOG_G_SUN = 3.7617, 4.4380

# per-feature scales for the medoid distance, from iso_compare_m15.py
_MEDOID_SCALES = {'Log_G': 0.06, 'Log_T': 0.02, 'Log_L': 0.06, 'Log_R': 0.06}

# physical bounds for the post-filter (overridable from the CLI)
LOGL_MIN, LOGL_MAX = -4.0, 3.6

ISO_COLUMNS = ["EEP","isochrone_age_yr","initial_mass","star_mass","star_mdot",
    "he_core_mass","c_core_mass","log_L","log_LH","log_LHe","log_Teff",
    "log_R","log_g","surface_h1","surface_he3","surface_he4","surface_c12",
    "surface_o16","log_center_T","log_center_Rho","center_gamma",
    "center_h1","center_he4","center_c12"]


def phys_from_varfile(path):
    """22-element physics block, read BY LINE POSITION.

    parse_varfeh() in iso_compare_m15.py keys off the second whitespace token,
    which breaks on our files because line 2 reads 'He Abundance' (two words)
    and yields the key 'He' rather than 'Yprim'.  Position is unambiguous.
    Element 17 (C12+He4 / deboercoef) is 1.0, matching physics_block_from_varfeh.
    """
    v = [float(l.split()[0]) for l in open(path) if l.strip()]
    if len(v) < 21:
        raise ValueError(f'{path}: expected >=21 lines, got {len(v)}')
    return np.array([
        v[0], v[1], v[2],            # FeH, Y, [a/Fe]
        v[3], v[4], v[5],            # mixing length, FGRY, FGRZ
        v[6], v[7], v[8],            # KTTAU, ALPHAE, ALPHAC
        v[9], v[10], v[11], v[12], v[13], v[14], v[15],   # SSTANDARD 1-7
        1.0,                         # C12+He4 (deboercoef)
        v[16], v[17], v[18], v[19], v[20],                # alex, opal2, talpha, plas, co
    ], dtype=np.float64)


def load_emulator(device):
    fmt = json.load(open(FMT))
    n_feat = sum(2 if s == 'mantissa_exponent' else 1 for s in fmt['data_scaler'])
    n_ctx  = sum(2 if s == 'mantissa_exponent' else 1 for s in fmt['label_scaler'])
    flow = zuko.flows.NSF(features=n_feat, context=n_ctx,
                          transforms=10, hidden_features=[256]*10).to(device)
    ck = torch.load(CKPT, map_location=device)
    # this checkpoint stores the state dict under 'model', with every key
    # prefixed 'flow.' from the training wrapper (and '_orig_mod.' if the run
    # used torch.compile, which config says it did)
    sd = ck.get('model', ck.get('model_state_dict', ck))
    clean = {}
    for k, v in sd.items():
        k = k.replace('_orig_mod.', '')
        if k.startswith('flow.'):
            k = k[len('flow.'):]
        clean[k] = v
    flow.load_state_dict(clean)
    flow.eval()
    return flow, fmt, n_feat


def _medoid(samples, headers, device):
    """Representative draw per label: the sample closest to the per-label median,
    distances scaled per feature. Same as representative_sample() upstream."""
    med = samples.median(dim=0).values
    diff = samples - med.unsqueeze(0)
    sc = torch.tensor([_MEDOID_SCALES.get(h, 1.0) for h in headers],
                      device=device, dtype=samples.dtype)
    d2 = ((diff / sc.view(1, 1, -1)) ** 2).sum(dim=2)
    idx = d2.argmin(dim=0)
    return samples[idx, torch.arange(samples.shape[1], device=device), :]


@torch.no_grad()
def emulate_central(flow, fmt, n_feat, phys, ages_yr, masses, device,
                    batch=200_000, n_draws=0):
    """n_draws=0 -> deterministic central point (z=0).
    n_draws>0 -> medoid of that many draws per label."""
    A, M = np.meshgrid(ages_yr, masses, indexing='ij')
    labels = np.hstack((np.tile(phys, (A.size, 1)),
                        A.reshape(-1, 1), M.reshape(-1, 1))).astype(np.float32)
    # sampling expands each label into n_draws rows, so the memory that matters
    # is batch * n_draws. Scale the label batch down so that product stays put;
    # at 64 draws an unscaled batch of 100k becomes 6.4M rows and OOMs a 24 GB card.
    eff = max(1000, batch // max(1, n_draws))
    batch = eff if n_draws > 0 else batch
    out = []
    for i in range(0, len(labels), batch):
        chunk = torch.tensor(labels[i:i+batch])
        norm = dt.norm_data_torch(chunk, fmt['label_scaler'], fmt['label_scale']).to(device)
        d = flow(norm)
        if n_draws > 0:
            raw = d.sample((n_draws,))               # (n_draws, B, n_feat)
            ns, nb, nf = raw.shape
            phys_s = dt.inverse_norm_data_torch(
                raw.reshape(ns * nb, nf), fmt['data_scaler'],
                fmt['data_scale']).reshape(ns, nb, nf)
            x = _medoid(phys_s, fmt['data_header'], device)
            out.append(x.cpu().numpy())
        else:
            z0 = torch.zeros(norm.shape[0], n_feat, device=device, dtype=norm.dtype)
            x = d.transform.inv(z0)                  # central, no sampling
            out.append(dt.inverse_norm_data_torch(
                x, fmt['data_scaler'], fmt['data_scale']).cpu().numpy())
    arr = np.vstack(out)
    df = pd.DataFrame(arr, columns=fmt['data_header'])
    df['isochrone_age_yr'] = A.reshape(-1)
    df['initial_mass'] = M.reshape(-1)
    # physical closure: R and g from L and Teff (same as iso_compare_m15.py)
    logR = 0.5*(df['Log_L'].values - 4.0*(df['Log_T'].values - LOG_T_SUN))
    df['Log_R'] = logR
    df['Log_G'] = LOG_G_SUN + np.log10(df['initial_mass'].values) - 2.0*logR
    # Same physical sanity cut as iso_compare_m15.py. The flow is poorly
    # constrained at the very lowest masses and emits a small number of absurd
    # points there (log L down to -18 around 0.15-0.20 Msun, ~1% of rows);
    # they sit far below the sCMD fit window but must not reach the files.
    # Upstream uses (-10, 10) in log L, which is far too loose: it lets through
    # log L = -18 at the lowest masses and +5.7 at masses that are already dead
    # at the isochrone's age (the fixed mass grid runs past the turnoff at old
    # ages). Real NGC 6791 isochrones span about -2.3 to +3.4, so bracket that
    # with margin instead.
    keep = ((df['Log_T'] > 3) & (df['Log_T'] < 4.5)
            & (df['Log_L'] > LOGL_MIN) & (df['Log_L'] < LOGL_MAX))
    n_drop = int((~keep).sum())
    if n_drop:
        print(f"    dropped {n_drop} unphysical rows ({n_drop/len(df)*100:.2f}%)")
    return df[keep].reset_index(drop=True)


def write_mist(df, path, feh, y, z, afe, n_ages):
    """MIST-format file the sCMD pipeline can read unchanged (11 '#' header lines)."""
    out = pd.DataFrame(0.0, index=np.arange(len(df)), columns=ISO_COLUMNS)
    out['EEP'] = np.arange(1, len(df)+1)
    out['isochrone_age_yr'] = df['isochrone_age_yr'].values
    out['initial_mass'] = df['initial_mass'].values
    out['star_mass'] = df['initial_mass'].values
    out['log_L'] = df['Log_L'].values
    out['log_Teff'] = df['Log_T'].values
    out['log_R'] = df['Log_R'].values
    out['log_g'] = df['Log_G'].values
    if 'M_He_Core' in df: out['he_core_mass'] = df['M_He_Core'].values
    if 'Y_Core' in df:    out['center_he4']   = df['Y_Core'].values
    npts = len(df)//max(n_ages,1)
    hdr = ("# MIST version number  = emul    \n"
           "# MESA revision number =    15140\n"
           "# " + "-"*86 + "\n"
           "#  Yinit        Zinit   [Fe/H]   [a/Fe]  v/vcrit\n"
           f"# {y:.4f}  {z:.5E}     {feh:.2f}     {afe:.2f}     0.00\n"
           "# " + "-"*86 + "\n"
           f"# number of isochrones = {n_ages:5d}\n"
           "# " + "-"*86 + "\n"
           f"# number of EEPs, cols = {npts:5d}   24\n"
           "# " + "  ".join(f"{i+1:>28d}" for i in range(24)) + "\n"
           "# " + "  ".join(f"{c:>28s}" for c in ISO_COLUMNS) + "\n")
    with open(path, 'w') as fh:
        fh.write(hdr)
        np.savetxt(fh, out.values, fmt='%28.16E')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('runs', nargs='+')
    ap.add_argument('--vardir', default='/dartfs-hpc/rc/home/8/f007058/mcdsep/MC_var_noco')
    ap.add_argument('--outdir', default='/dartfs-hpc/scratch/f007058/Isochrones_emul')
    ap.add_argument('--device', default='cuda' if torch.cuda.is_available() else 'cpu')
    ap.add_argument('--mass-low', type=float, default=0.15)
    ap.add_argument('--mass-high', type=float, default=1.60)
    ap.add_argument('--n-mass', type=int, default=1200)
    ap.add_argument('--n-draws', type=int, default=0,
                    help='0 = central (z=0, deterministic); >0 = medoid of N draws')
    ap.add_argument('--batch', type=int, default=200_000)
    ap.add_argument('--seed', type=int, default=0)
    ap.add_argument('--logl-min', type=float, default=-4.0)
    ap.add_argument('--logl-max', type=float, default=3.6)
    args = ap.parse_args()

    global LOGL_MIN, LOGL_MAX
    LOGL_MIN, LOGL_MAX = args.logl_min, args.logl_max
    ages = np.round(np.arange(4.0, 12.0001, 0.1), 4) * 1e9
    masses = np.linspace(args.mass_low, args.mass_high, args.n_mass)
    flow, fmt, n_feat = load_emulator(args.device)
    mode = f'medoid of {args.n_draws} draws' if args.n_draws > 0 else 'central (z=0)'
    print(f"emulator loaded on {args.device}; {len(ages)} ages x {len(masses)} masses; {mode}")

    import subprocess
    for run in args.runs:
        vf = f'{args.vardir}/varfeh+031.{run}'
        if not os.path.exists(vf):
            print(f"  {run}: no var file"); continue
        phys = phys_from_varfile(vf)
        feh, y, afe = phys[0], phys[1], phys[2]
        zz = float(subprocess.run(['/dartfs-hpc/rc/home/8/f007058/mcdsep/run/gs98getz',
                                   str(feh), str(afe)], capture_output=True,
                                  text=True).stdout.strip() or 0.0)
        torch.manual_seed(args.seed)
        df = emulate_central(flow, fmt, n_feat, phys, ages, masses, args.device,
                             batch=args.batch, n_draws=args.n_draws)
        d = f'{args.outdir}/{run}/isos'; os.makedirs(d, exist_ok=True)
        p = f'{d}/isochrones_example.txt'
        write_mist(df, p, feh, y, zz, afe, len(ages))
        print(f"  {run} -> {os.path.getsize(p)/1e6:.1f} MB  "
              f"logL [{df['Log_L'].min():+.2f},{df['Log_L'].max():+.2f}]")

if __name__ == '__main__':
    main()
