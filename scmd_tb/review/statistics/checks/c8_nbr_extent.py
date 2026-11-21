import sys, numpy as np, pandas as pd
sys.path.insert(0, "/dartfs-hpc/scratch/f007058/scmd_tb/review/statistics/dev")
from isochronetoolbox.scmd.ast import NearestErrors, centre
A = pd.read_csv("/dartfs-hpc/scratch/f007058/Files_from_Martin/47Tuc_data/47Tucartstars.dat", dtype=np.float32)
zpv, zpi = 30.5742, 29.6777
x, y = A.InputX.values, A.InputY.values; xc, yc = centre(x, y); r = np.hypot(x - xc, y - yc)
vin = A.InputF606W.values.astype(float); iin = A.InputF814W.values.astype(float)
win = (vin + zpv > 15.02) & (vin + zpv < 19.02)
for metric, sp in (("linear", None), ("crowding", (16.48 - zpv, 15.55 - zpi))):
    ne = NearestErrors(A, k=30, metric=metric, split=sp)
    nb = ne.nbr[win]
    dv = np.abs(vin[nb] - vin[win][:, None]).max(1); dr = np.abs(r[nb] - r[win][:, None]).max(1)
    out = []
    for lo, hi in [(0, 200), (200, 600), (600, 1500), (1500, 5000)]:
        s = (r[win] >= lo) & (r[win] < hi)
        out.append(f"r{lo}-{hi}: max|dV| med {np.median(dv[s]):.3f} max|dr| med {np.median(dr[s]):.0f}px")
    print(f"{metric:8s} k=30 neighbourhood extent: " + " | ".join(out))
