#!/usr/bin/env python3
"""P22 Amendment 12 analysis (Gaspar-Kukal arms). Usage: analyze_amend12.py [--smoke]. Reads results/raw_amend12_{b8,b32,b128}.csv, cells_meta_amend12_*.csv;
writes results/amend12_table.csv, amend12_minL.csv, amend12_report.txt. Bars are those of Amendment 9 / PLAN Amendment 12; verdict lines printed literally."""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import beta
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
SFX = "_smoke" if "--smoke" in sys.argv else ""
G = {"aniso3_w8": [32, 48, 64, 80, 96, 112, 128, 160, 192, 256], "aniso2_w16": [128, 160, 192, 224, 256, 320, 384, 512, 768]}
REAL = lambda fam, l: (fam == "aniso2_w16" and l >= 128) or (fam == "aniso3_w8" and l >= 48)
GLOBAL = ("GKW", "GKL", "GKWm", "GKWb", "E0g"); TWO = ("GK2W", "GK2L", "GK2Wb", "E0")
LINES = []
def P(s=""): print(s); LINES.append(str(s))
def cp(k, n):
    k, n = int(k), int(n)
    if n == 0: return (np.nan, np.nan)
    return (beta.ppf(.025, k, n - k + 1) if k > 0 else 0.0, beta.ppf(.975, k + 1, n - k) if k < n else 1.0)
def fmt(k, n):
    lo, hi = cp(k, n); return "%d/%d = %.1f%% (CP %.1f-%.1f)" % (k, n, 100 * k / max(n, 1), 100 * lo, 100 * hi)

def load(g):
    p = os.path.join(RES, f"raw_amend12_{g}{SFX}.csv"); m = os.path.join(RES, f"cells_meta_amend12_{g}{SFX}.csv")
    if not (os.path.exists(p) and os.path.exists(m)): return None
    d = pd.read_csv(p, low_memory=False); meta = pd.read_csv(m)
    d = d.merge(meta[["id", "fam", "L", "kind"]], left_on="cell", right_on="id", how="inner")
    for c in ("decl_neg", "decl_pos", "cov_nom"):
        if c not in d.columns: d[c] = np.nan
    d["decl_any"] = ((d.decl_neg.fillna(0) + d.decl_pos.fillna(0)) > 0).astype(int)
    d["real"] = [REAL(f, l) and k == "crossover" for f, l, k in zip(d.fam, d.L, d.kind)]
    d["valid"] = d.abstain == 0
    return d

def calib(x, two):
    c = x[x.kind == "control"]; ok = c[c.valid]
    k, n = int(ok.cov_simult.sum()), len(ok)
    per = ok.groupby("cell").cov_simult.agg(["sum", "count"]); badc = [int(i) for i, r in per.iterrows() if cp(r["sum"], r["count"])[1] < .93]
    C = n > 0 and k / n >= .93 and not badc
    kn = int(ok.cov_nom.sum()); r = dict(est=x.est.iloc[0], n=n, abst=round(1 - n / max(len(c), 1), 3), C=bool(C), cov=fmt(k, n), cell_cp_upper_lt_93=badc, cov_nominal=(fmt(kn, n) if x.est.iloc[0].startswith("GK") else "NA"),
                                         bias=round(float(ok.bias.mean()), 4), hw=round(float(ok.hw.mean()), 4))
    rr = x[x.real & x.valid]; r["cov_realised_crossover(X-descr)"] = fmt(int(rr.cov_simult.sum()), len(rr))
    if two:
        fk = int(ok.decl_any.sum()); pf = ok.groupby("cell").decl_any.agg(["sum", "count"]); badf = [int(i) for i, q in pf.iterrows() if cp(q["sum"], q["count"])[0] > .05]
        r.update(F=bool(n > 0 and fk / n <= .05 and not badf), FPR=fmt(fk, n), F_bad=badf)
    return r

def minL(x):
    cx = x[x.kind == "crossover"].copy(); cx["det"] = ((cx.valid) & (cx.decl_neg == 1)).astype(float)
    pw = cx.groupby(["fam", "L"]).det.agg(["sum", "count"]).reset_index(); pw["power"] = pw["sum"] / pw["count"]; out = {}
    for fam, grid in G.items():
        g = pw[pw.fam == fam]; Pw = dict(zip(g.L, g.power)); idx = len(grid)
        for i in range(len(grid)):
            if all((l in Pw) and Pw[l] >= .8 for l in grid[i:]): idx = i; break
        out[fam] = grid[idx] if idx < len(grid) else np.inf
        out[fam + "_power"] = " ".join(f"{int(l)}:{p:.2f}" for l, p in zip(g.L, g.power))
    return out

rows = []; mrows = []
for g, w in (("b32", 32), ("b8", 8), ("b128", 128)):
    d = load(g)
    if d is None: P(f"[wpn {w}] no data"); continue
    P(f"\n=== wpn {w}: {d.cell.nunique()} cells, {len(d)} rows ===")
    for est in GLOBAL + TWO:
        x = d[d.est == est]
        if len(x) == 0: continue
        r = calib(x, est in TWO); r["wpn"] = w; rows.append(r); P(" ".join(f"{k}={v}" for k, v in r.items()))
        if est in TWO:
            m = minL(x); m.update(est=est, wpn=w); mrows.append(m); P(f"   minL: aniso2_w16={m['aniso2_w16']} aniso3_w8={m['aniso3_w8']}   power2D {m['aniso2_w16_power']} | 3D {m['aniso3_w8_power']}")
pd.DataFrame(rows).to_csv(os.path.join(RES, f"amend12_table{SFX}.csv"), index=False); pd.DataFrame(mrows).to_csv(os.path.join(RES, f"amend12_minL{SFX}.csv"), index=False)
P("\nVERDICT LINES (wpn 32, gated)")
t = {(r["est"], r["wpn"]): r for r in rows}
for e in GLOBAL:
    if (e, 32) in t: P(f"{e}: C {'PASS' if t[(e,32)]['C'] else 'FAIL'}  ({t[(e,32)]['cov']})")
for e in TWO:
    if (e, 32) in t: P(f"{e}: C {'PASS' if t[(e,32)]['C'] else 'FAIL'}, F {'PASS' if t[(e,32)]['F'] else 'FAIL'}  ({t[(e,32)]['cov']}; FPR {t[(e,32)]['FPR']})")
open(os.path.join(RES, f"amend12_report{SFX}.txt"), "w").write("\n".join(LINES) + "\n")
