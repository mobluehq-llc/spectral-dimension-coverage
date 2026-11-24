#!/usr/bin/env python3
"""ROUND 5 (Gemini r4#? / ChatGPT r5 s3): (a) L_min-bracketing counts with CP intervals for the Table II arms E1q, E4c, E4r, E6, Mb, CL
(mirrors analyze_amend4.py / analyze_amend5.py: power = decl_neg over ALL 200 datasets of each realised crossover cell, abstention a miss);
(b) G3/G4 gate inputs of the holdout (from results_amend10_11/tables_11_cells.csv). Reads stored per-dataset results only; nothing is re-run.
Usage: round5_tables.py --out DIR      Fails loudly (MISMATCH) if a recomputed L_min differs from Table II."""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import beta
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results"); R11 = os.path.join(HERE, "results_amend10_11")
OUT = sys.argv[sys.argv.index("--out") + 1]; os.makedirs(OUT, exist_ok=True)
def cp(k, n):
    k, n = int(k), int(n)
    return (beta.ppf(.025, k, n - k + 1) if k > 0 else 0.0, beta.ppf(.975, k + 1, n - k) if k < n else 1.0)
meta = pd.read_csv(os.path.join(RES, "cells_meta.csv"), keep_default_na=False); meta["note"] = meta["note"].astype(str)
meta = meta[(meta.note == "") & meta.kind.isin(["control", "crossover"])][["id", "fam", "L", "N", "kind", "max_drop"]]
m5 = pd.read_csv(os.path.join(RES, "cells_meta_amend5.csv"))[["id", "fam", "L", "N", "kind", "max_drop"]]
mb = pd.read_csv(os.path.join(RES, "cells_meta_amend6_big2d.csv"))[["id", "fam", "L", "N", "kind", "max_drop"]]  # 2D sizes 384/512/768 for E4c, E4r, E6
meta = pd.concat([meta, m5[m5.id >= 100], mb], ignore_index=True); meta["max_drop"] = pd.to_numeric(meta.max_drop, errors="coerce")
SRC = {"E1q": "raw_amend.csv", "E4c": "raw_amend3.csv", "E4r": "raw_amend3.csv", "E6": "raw_amend4.csv", "Mb": "raw_amend5.csv", "CL": "raw_amend5.csv"}
# Table II values: (2D L_min, 3D L_min); None = not reached (E4c 2D, sweep ended at 768)
EXP = {"E1q": (192, 96), "E4c": (None, 192), "E4r": (768, 192), "E6": (256, 96), "Mb": (256, 128), "CL": (128, 64)}
cache = {}; rows = []; bad = []
for est, f in SRC.items():
    if f not in cache: cache[f] = pd.read_csv(os.path.join(RES, f), low_memory=False)
    d = cache[f]
    if est in ("E4c", "E4r", "E6"): d = pd.concat([d, pd.read_csv(os.path.join(RES, "raw_amend6_big2d.csv"), low_memory=False)], ignore_index=True)
    d = d[d.est == est].merge(meta, left_on="cell", right_on="id")
    d = d[(d.kind == "crossover") & (d.max_drop >= 0.6)]
    for fam, ex in zip(("aniso2_w16", "aniso3_w8"), EXP[est]):
        t = d[d.fam == fam].groupby("L").apply(lambda x: pd.Series(dict(k=x.decl_neg.fillna(0).sum(), n=len(x))), include_groups=False).sort_index()
        pw = (t.k / t.n).values
        ip = next((i for i in range(len(t)) if (pw[i:] >= 0.8).all()), None)
        got = int(t.index[ip]) if ip is not None else None
        if got != ex: bad.append(f"{est} {fam}: want {ex} got {got}")
        ref = ip if ip is not None else len(t) - 1
        for i in range(max(ref - 1, 0), min(ref + 2, len(t))) if ip is not None else range(max(ref - 2, 0), ref + 1):
            k, n = int(t.k.iloc[i]), int(t.n.iloc[i]); lo, hi = cp(k, n)
            rows.append(dict(arm=est, fam=fam, L=int(t.index[i]), k=k, n=n, pct=100 * k / n, lo=100 * lo, hi=100 * hi, is_Lmin=(ip is not None and i == ip), Lmin_reported=ex))
B = pd.DataFrame(rows); B.to_csv(os.path.join(OUT, "r5_Lmin_brackets.csv"), index=False)
if bad: print("MISMATCH", bad); sys.exit(1)
print("L_min recomputed equals Table II for all 12 arm-family entries")
with open(os.path.join(OUT, "tab_r5_brackets.tex"), "w") as fh:
    for (arm, fam), g in B.groupby(["arm", "fam"], sort=False):
        fl = "2D" if fam == "aniso2_w16" else "3D"
        for j, r in enumerate(g.itertuples()):
            mk = "$\\leftarrow L_{\\min}$" if r.is_Lmin else ""
            fh.write(f"{arm if j == 0 else ''} & {fl if j == 0 else ''} & {r.L} & {r.k}/{r.n} & {r.pct:.1f} ({r.lo:.1f}--{r.hi:.1f}) & {mk} \\\\\n")
        fh.write("\\hline\n")
# (b) holdout gate inputs
a = pd.read_csv(os.path.join(R11, "tables_11_cells.csv")); G = []
cx = a[(a.kind == "crossover") & a.real.astype(bool)]; ct = a[a.kind == "control"]
for est in ("E4c", "E4r"):
    for w in (8, 32, 128):
        y = cx[(cx.wpn == w)]; k = int(y[f"{est}_cov_k"].sum()); n = int(y[f"{est}_nvalid"].sum()); lo, hi = cp(k, n)
        G.append(dict(gate="G3", arm=est, wpn=w, k=k, n=n, nn=int(len(y) * 200), pct=100 * k / n, lo=100 * lo, hi=100 * hi, worst_cell="", worst_hi=np.nan))
for est in ("E4c", "E4r", "E7", "E8"):
    for w in (8, 32, 128):
        y = ct[(ct.wpn == w)]; k = int(y[f"{est}_cov_k"].sum()); n = int(y[f"{est}_nvalid"].sum()); lo, hi = cp(k, n)
        his = [(cp(r[f"{est}_cov_k"], r[f"{est}_nvalid"])[1], f"{r['fam']} {int(r['L'])}") for _, r in y.iterrows() if r[f"{est}_nvalid"] > 0]
        wh = min(his)
        G.append(dict(gate="G4", arm=est, wpn=w, k=k, n=n, nn=int(len(y) * 200), pct=100 * k / n, lo=100 * lo, hi=100 * hi, worst_cell=wh[1], worst_hi=100 * wh[0]))
G = pd.DataFrame(G); G.to_csv(os.path.join(OUT, "r5_gate_inputs.csv"), index=False)
print(G.to_string())

with open(os.path.join(OUT, "tab_r5_gates.tex"), "w") as fh:
    for r in G.itertuples():
        w = "" if r.gate == "G3" else f"{r.worst_cell} ({r.worst_hi:.1f})"
        fh.write(f"{r.gate} & {r.arm} & {r.wpn} & {r.k}/{r.n} & {r.pct:.1f} ({r.lo:.1f}--{r.hi:.1f}) & {w} \\\\\n")
        if r.wpn == 128: fh.write("\\hline\n")
