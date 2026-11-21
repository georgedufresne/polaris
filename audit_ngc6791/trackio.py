"""Readers for DSEP MESA-format track files and MIST-format isochrones."""
import gzip, io, numpy as np

TRACK_COLS = ['star_age','star_mass','star_mdot','he_core_mass','c_core_mass',
              'log_L','log_LH','log_LHe','log_Teff','log_R','log_g',
              'surface_h1','surface_he3','surface_he4','surface_c12','surface_o16',
              'log_center_T','log_center_Rho','center_gamma',
              'center_h1','center_he4','center_c12']

def read_track(path):
    """DSEP writes MESA history with ONE VALUE PER LINE. Returns dict of arrays."""
    op = gzip.open if str(path).endswith('.gz') else open
    with op(path, 'rt', errors='replace') as fh:
        lines = fh.read().split('\n')
    # locate the 22-name column header, data follows
    start = None
    for i, l in enumerate(lines):
        if 'star_age' in l and 'log_Teff' in l:
            start = i + 1; break
    if start is None:
        return None
    vals = []
    for l in lines[start:]:
        s = l.strip()
        if not s:
            continue
        try:
            vals.append(float(s))
        except ValueError:
            continue
    n = len(vals) // len(TRACK_COLS)
    if n == 0:
        return None
    a = np.array(vals[:n*len(TRACK_COLS)]).reshape(n, len(TRACK_COLS))
    return {c: a[:, j] for j, c in enumerate(TRACK_COLS)}

def read_mist_iso(path):
    """MIST-format isochrone file -> (colnames, array, header dict)."""
    with open(path) as fh:
        lines = fh.read().split('\n')
    cols = None; rows = []; hdr = {}
    for l in lines:
        if l.startswith('#'):
            if 'EEP' in l and 'isochrone_age_yr' in l:
                cols = l.lstrip('#').split()
            continue
        s = l.split()
        if not s or cols is None:
            continue
        try:
            rows.append([float(x) for x in s])
        except ValueError:
            continue
    if not rows:
        return None, None, hdr
    return cols, np.array(rows), hdr

def read_dsep_iso(path):
    """DSEP native isochrone file -> {age_myr: dict of arrays}.

    Layout, repeated per age block:
      #NPTS MIX-LEN  OVERSH AGE(MYR)   Y   Z   ZEFF  [Fe/H] [alpha/Fe]
      #355 1.088800  0.0000  6000.  0.2516 ...
      #EEP  MASS  Log G  Log Teff  Log(L/Ls)  NUM0 NUM-2.35 NUM-4.0  V  V-I  F606W  F606W-F814W
      284* 0.709811  4.567819  3.698668 -0.5365753 ...
    """
    import numpy as np
    out={}; age=None; rows=[]
    for line in open(path):
        if line.startswith('#'):
            p=line[1:].split()
            if p and p[0].replace('.','',1).isdigit() and len(p)>=4:
                if age is not None and rows: out[age]=rows
                try: age=float(p[3].rstrip('.')); rows=[]
                except (ValueError,IndexError): pass
            continue
        p=line.split()
        if len(p)<5: continue
        try: rows.append([float(p[0].rstrip('*')),float(p[1]),float(p[2]),float(p[3]),float(p[4])])
        except ValueError: continue
    if age is not None and rows: out[age]=rows
    return {a:{'EEP':np.array([r[0] for r in v]),
               'mass':np.array([r[1] for r in v]),
               'Log_G':np.array([r[2] for r in v]),
               'Log_T':np.array([r[3] for r in v]),
               'Log_L':np.array([r[4] for r in v])} for a,v in out.items()}
