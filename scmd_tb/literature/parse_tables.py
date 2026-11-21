#!/usr/bin/env python3
"""Raw literature tables (as pasted) -> tidy CSVs."""
import re
import pandas as pd

V = r"(—|-?\d+\.\d+±\d+\.\d+)"
def pm(s):
    return (None, None) if s == "—" else tuple(float(x) for x in s.split("±"))

# Milone et al. 2012, Table 2
rows, cl = [], None
pat = re.compile(rf"^(?P<head>.*?)(?P<region>rC sample|rC−HM sample|roHM sample|WFC field) {V} {V} {V} {V}$")
for line in open("milone2012_table2_raw.txt", encoding="utf-8"):
    m = pat.match(line.strip())
    if not m:
        continue
    head, region = m["head"].strip(), m["region"]
    rmin = None
    if region == "WFC field":
        rmin = float(re.match(r"Rmin=(\d+\.\d+)", head).group(1))
    elif region == "rC sample":
        cl = head
    r = {"cluster": cl, "region": region, "rmin": rmin}
    for k, s in zip(["f_q05", "f_q06", "f_q07", "f_tot"], m.groups()[2:]):
        r[k], r[k + "_err"] = pm(s)
    rows.append(r)
b = pd.DataFrame(rows)
assert len(b) == 4 * b.cluster.nunique(), (len(b), b.cluster.nunique())
b.to_csv("milone2012_binary_fractions.csv", index=False)
print(f"milone2012_binary_fractions.csv: {b.cluster.nunique()} clusters, {len(b)} rows")

# mass-function table
cols = ["alpha", "log_mlum", "log_mdyn", "rh_pc", "f_remn", "trh_gyr", "mdyn_lv", "rj_pc"]
out = []
for line in open("gc_mass_function_table_raw.txt", encoding="utf-8"):
    m = re.match(r"^(\d{3,4})\s?([−-]?\d.*)$", line.strip())
    if not m:
        continue
    toks = m[2].replace("−", "-").replace(" ± ", "±").split()
    r = {"ngc": int(m[1])}
    if len(toks) == 3:                                   # NGC 6584: alpha, log Mlum, rh only
        keys = ["alpha", "log_mlum", "rh_pc"]
    else:
        keys = cols
    for k, t in zip(keys, toks):
        if "±" in t:
            r[k], r[k + "_err"] = map(float, t.split("±"))
        else:
            r[k] = float(t)
    out.append(r)
order = ["ngc"] + [c for k in cols for c in (k, k + "_err") if c not in ("trh_gyr_err", "rj_pc_err")]
mf = pd.DataFrame(out).reindex(columns=order)
mf.to_csv("gc_mass_functions.csv", index=False)
print(f"gc_mass_functions.csv: {len(mf)} clusters")
print(b[b.cluster.isin(["NGC 104", "NGC 7078"])].to_string(index=False))
print(mf[mf.ngc.isin([104, 7078])].to_string(index=False))
