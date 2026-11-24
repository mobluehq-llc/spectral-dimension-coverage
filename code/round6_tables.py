#!/usr/bin/env python3
"""ROUND 6 (ChatGPT r6 s3): bracketing counts and CP intervals for every retained threshold of Table II not already in
Tables S10/S11 (E1, E1x, E1u, E3, M) and for every entry of Table IV (E0, E1q, E6, E4c, E4r at 8/32/128 walkers per node).
Power = declared drops over ALL 200 datasets of each realised crossover cell, abstention a miss. Reads stored per-dataset
results only; nothing is re-run. Usage: round6_tables.py --out DIR. Fails loudly (MISMATCH) if a recomputed L_min differs
from the reported one. Also writes sweep_limits.csv (largest tested L per arm and family)."""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import beta
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
OUT = sys.argv[sys.argv.index("--out") + 1]; os.makedirs(OUT, exist_ok=True)
def cp(k, n):
    k, n = int(k), int(n)
    return (beta.ppf(.025, k, n - k + 1) if k > 0 else 0.0, beta.ppf(.975, k + 1, n - k) if k < n else 1.0)
FAM = ("aniso2_w16", "aniso3_w8"); bad = []; rows = []; lim = []
def lmin(pw):
    return next((i for i in range(len(pw)) if (pw[i:] >= 0.8).all()), None)
def window(ip, n):
    if ip is not None: return list(range(max(ip - 1, 0), min(ip + 2, n)))
    return list(range(max(n - 3, 0), n))
# ---- Table II arms (32 walkers per node)
m_old = pd.read_csv(os.path.join(RES, "cells_meta.csv"), keep_default_na=False); m_old["note"] = m_old["note"].astype(str)
m_old = m_old[(m_old.note == "") & m_old.kind.isin(["control", "crossover"])][["id", "fam", "L", "N", "kind", "max_drop"]]
m5 = pd.read_csv(os.path.join(RES, "cells_meta_amend5.csv"))[["id", "fam", "L", "N", "kind", "max_drop"]]
meta = pd.concat([m_old, m5[m5.id >= 100]], ignore_index=True); meta["max_drop"] = pd.to_numeric(meta.max_drop, errors="coerce")
raw = pd.concat([pd.read_csv(os.path.join(RES, "raw.csv")), pd.read_csv(os.path.join(RES, "raw_amend.csv")),
                 pd.read_csv(os.path.join(RES, "raw_amend5.csv"), low_memory=False)], ignore_index=True)
EXP = {"E1": (256, 96), "E1x": (192, 96), "E1u": (None, 96), "E3": (128, 48), "M": (192, 128)}
for est, ex in EXP.items():
    d = raw[raw.est == est].merge(meta, left_on="cell", right_on="id")
    d = d[(d.kind == "crossover") & (d.max_drop >= 0.6)]
    for fam, e in zip(FAM, ex):
        t = d[d.fam == fam].groupby("L").apply(lambda x: pd.Series(dict(k=x.decl_neg.fillna(0).sum(), n=len(x))), include_groups=False).sort_index()
        ip = lmin((t.k / t.n).values); got = int(t.index[ip]) if ip is not None else None
        if got != e: bad.append(f"{est} {fam}: want {e} got {got}")
        lim.append(dict(arm=est, fam=fam, budget=32, max_L=int(t.index.max())))
        rows.append(dict(tab="II", arm=est, fam=fam, wpn=32, Lmin=got, pts=[(int(t.index[i]), int(t.k.iloc[i]), int(t.n.iloc[i])) for i in window(ip, len(t))]))
# ---- Table IV (all arms, budgets)
P = pd.read_csv(os.path.join(RES, "amend7_power.csv")); MM = pd.read_csv(os.path.join(RES, "amend7_minL.csv"))
for (est, fam, w), g in P.groupby(["est", "fam", "wpn"], sort=False):
    g = g.sort_values("L"); k = (g.power * g.n).round().astype(int).values
    assert np.allclose(k, (g.power * g.n).values, atol=1e-6), "power*n not an integer"
    ip = lmin((k / g.n.values)); got = int(g.L.iloc[ip]) if ip is not None else None
    r = MM[(MM.est == est) & (MM.fam == fam) & (MM.wpn == w)].iloc[0]
    want = int(r.meas_minL) if str(r.meas_minL).strip().isdigit() else None  # '>768' etc. = not reached
    if got != want: bad.append(f"IV {est} {fam} {w}: want {want} got {got}")
    lim.append(dict(arm=est, fam=fam, budget=w, max_L=int(g.L.max())))
    rows.append(dict(tab="IV", arm=est, fam=fam, wpn=w, Lmin=got, pts=[(int(g.L.iloc[i]), int(k[i]), int(g.n.iloc[i])) for i in window(ip, len(g))]))
pd.DataFrame(lim).to_csv(os.path.join(OUT, "sweep_limits.csv"), index=False)
R = pd.DataFrame(rows); R.to_csv(os.path.join(OUT, "r6_brackets.csv"), index=False)
if bad: print("MISMATCH", bad); sys.exit(1)
print("L_min recomputed equals the reported value for all", len(R), "arm-family-budget entries")
def cell(L, k, n, mark):
    lo, hi = cp(k, n); return f"{L}: {k}/{n} ({100*k/n:.1f}; {100*lo:.1f}--{100*hi:.1f})" + ("$^\\ast$" if mark else "")
def dump(tab, fn):
    with open(os.path.join(OUT, fn), "w") as fh:
        prev = None
        for r in R[R.tab == tab].itertuples():
            fl = "2D" if r.fam == "aniso2_w16" else "3D"
            cs = [cell(L, k, n, L == r.Lmin) for (L, k, n) in r.pts]; cs += [""] * (3 - len(cs))
            lm = int(r.Lmin) if not pd.isna(r.Lmin) else "$>%d$" % r.pts[-1][0]
            fh.write(f"{r.arm} & {fl} & {r.wpn} & {lm} & {cs[0]} & {cs[1]} & {cs[2]} \\\\\n")
dump("II", "tab_r6_brk2.tex"); dump("IV", "tab_r6_brk4.tex")
