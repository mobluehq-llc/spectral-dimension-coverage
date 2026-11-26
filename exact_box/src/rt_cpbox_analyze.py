#!/usr/bin/env python3
"""BD analysis: reads results/bd_*.jsonl (+ stored nodrop_data.csv for the nominal join) and writes results/BD_tables.md. Usage: rt_cpbox_analyze.py [results_dir]"""
import sys, os, json
import numpy as np, pandas as pd
from scipy.stats import beta
RES = sys.argv[1] if len(sys.argv) > 1 else "results"
def cpi(k, n):
    if n == 0: return (np.nan, np.nan)
    return (0.0 if k == 0 else beta.ppf(0.025, k, n - k + 1), 1.0 if k == n else beta.ppf(0.975, k + 1, n - k))
def fmt(k, n):
    lo, hi = cpi(k, n); return "%d/%d [%.3f, %.3f]" % (k, n, lo, hi)
def load(m):
    p = os.path.join(RES, "bd_%s.jsonl" % m)
    return pd.DataFrame([json.loads(l) for l in open(p)]) if os.path.exists(p) else pd.DataFrame()
out = []
def sec(t): out.append("\n## " + t + "\n")
h = os.path.join(RES, "bd_harness_summary.json")
if os.path.exists(h):
    s = json.load(open(h)); sec("Bar 1 harness (ring 1024, 32/node)"); out.append("covered %d/%d = %.4f, one-sided CP lower limit %.4f (needs >= 0.94): **%s**" % (s["covered"], s["n"], s["rate"], s["one_sided_cp_lower"], "PASS" if s["pass"] else "FAIL"))
errs = 0
a = load("a"); b = load("b"); c = load("c"); d = load("d"); e = load("e"); f = load("f")
for df in (a, b, c, d, e, f):
    if len(df) and "g_outcome" in df: errs += int((df.g_outcome == "ERROR").sum())
out.append("\nERROR rows across all modes: %d (must be 0 for the report to stand)" % errs)
if len(a):
    sec("Arm (a): same-seed rerun, CP box, declarations out of datasets (95% CP)")
    nom = pd.read_csv(os.path.join(os.path.dirname(RES.rstrip("/")) or ".", "results", "nodrop_data.csv")) if False else None
    rows = []
    for (fam, Lv, w), g in a.groupby(["fam", "L", "wpn"]):
        n = len(g); rows.append(dict(cell="%s L=%d" % (fam, Lv), wpn=w, n=n, true_drop=round(g.true_drop.iloc[0], 4),
            global_DROP=fmt(int((g.g_outcome == "DROP").sum()), n), joint_DROP=fmt(int((g.j_drop == 1).sum()), n) if "j_drop" in g else "",
            pointwise_DROP=fmt(int((g.e7_neg == 1).sum()), n) if "e7_neg" in g else "", secant_DROP=fmt(int((g.s_neg == 1).sum()), n) if "s_neg" in g else "",
            box_cover=fmt(int(g.box_cov.sum()), n), hw_ratio=round(g.hw_ratio.mean(), 3), sec_per_ds=round(g.wall.mean(), 1)))
    out.append(pd.DataFrame(rows).to_markdown(index=False))
    sp = os.path.join(RES, "nodrop_data.csv")
    if os.path.exists(sp):
        st = pd.read_csv(sp); m = a[a.wpn == 32].merge(st[["cell", "ds", "g_outcome", "j_drop"]], on=["cell", "ds"], how="inner", suffixes=("", "_nom"))
        out.append("\nJoin with the stored nominal-box outcome (32/node, same seeds):\n")
        t = []
        for cell, g in m.groupby("cell"):
            t.append(dict(cell=cell, fam=g.fam.iloc[0], n=len(g), nominal_global_DROP=int((g.g_outcome_nom == "DROP").sum()), cp_global_DROP=int((g.g_outcome == "DROP").sum()),
                          both=int(((g.g_outcome_nom == "DROP") & (g.g_outcome == "DROP")).sum()), cp_only=int(((g.g_outcome_nom != "DROP") & (g.g_outcome == "DROP")).sum())))
        out.append(pd.DataFrame(t).to_markdown(index=False))
if len(a) or len(b):
    sec("Bar 2: level (false declarations = DROP on a graph with exact drop < 0.245), CP box")
    fa = a[(a.true_drop < 0.245) & (a.g_outcome.isin(["DROP", "WITNESS", "INCONCLUSIVE", "ERROR", "BOXINFEASIBLE", "NOWINDOW"]))] if len(a) else pd.DataFrame()
    fb = b[(b.true_drop < 0.245) & (b.box == "cp")] if len(b) else pd.DataFrame()
    allf = pd.concat([fa, fb]); k = int((allf.g_outcome == "DROP").sum()); n = len(allf)
    lo, hi = cpi(k, n)
    out.append("pooled controls (a) + sub-margin graphs (b): false DROP %s; CP lower limit %.4f vs 0.05: **%s**" % (fmt(k, n), lo, "PASS (rate <= 5%)" if k / max(n, 1) <= 0.05 else ("CONTRADICTS section 2 (report as fault)" if lo > 0.05 else "observed above 5% but CP lower limit below 5%")))
    if len(a):
        out.append("\nBox coverage of exact q over all grid times, all (a) datasets: %s" % fmt(int(a.box_cov.sum()), len(a)))
if len(b):
    sec("Arm (b): near-null graphs, global DROP rate by box")
    rows = []
    for (cell, w, bx), g in b.groupby(["cell", "wpn", "box"]):
        rows.append(dict(cell=cell, L=g.L.iloc[0], N=g.N.iloc[0], exact_drop=round(g.true_drop.iloc[0], 4), wpn=w, box=bx, n=len(g), DROP=fmt(int((g.g_outcome == "DROP").sum()), len(g)), nowindow=int((g.g_outcome == "NOWINDOW").sum())))
    out.append(pd.DataFrame(rows).sort_values(["wpn", "box", "exact_drop"]).to_markdown(index=False))
if len(c) or len(a):
    sec("Arm (c): power map, global test, CP box (cells of (a) contribute their rows)")
    pm = pd.concat([x for x in (c, a[a.fam.isin(["aniso2_w16", "aniso3_w8"])] if len(a) else None) if x is not None and len(x)])
    pm = pm.drop_duplicates(["cell", "wpn", "ds"]); rows = []
    for (fam, Lv, w), g in pm.groupby(["fam", "L", "wpn"]):
        rows.append(dict(fam=fam, L=Lv, wpn=w, n=len(g), DROP=int((g.g_outcome == "DROP").sum()), CI=fmt(int((g.g_outcome == "DROP").sum()), len(g)), sec_per_ds=round(g.wall.mean(), 1)))
    tab = pd.DataFrame(rows).sort_values(["fam", "L", "wpn"]); out.append(tab.to_markdown(index=False))
    mono = []
    for (fam, Lv), g in tab.groupby(["fam", "L"]):
        g = g.sort_values("wpn"); ks = g.DROP.tolist(); mono.append("%s L=%d: %s %s" % (fam, Lv, ks, "monotone" if all(x <= y for x, y in zip(ks, ks[1:])) else "NOT monotone (check CIs)"))
    out.append("\nBar 4 monotonicity (raw counts at 8, 32, 128, 512):\n\n" + "\n".join("- " + x for x in mono))
    big = tab[(tab.fam == "aniso2_w16") & (tab.L == tab[tab.fam == "aniso2_w16"].L.max())]
    out.append("\nBar 3 (GO): C768xC16 DROP counts by wpn: %s" % dict(zip(big.wpn, big.DROP)))
if len(d):
    sec("Arm (d): zero-noise limit, exact q (1 +- eps), global test")
    out.append(d[["fam", "L", "wpn", "eps", "Jw", "g_outcome", "g_r", "g_mindelta", "g_nlp", "g_gate"]].sort_values(["fam", "L", "eps"]).to_markdown(index=False))
if len(e):
    sec("Arm (e): continuum certificates of CP-box DROPs")
    out.append(e.groupby(["cell", "wpn"]).verdict.value_counts().unstack(fill_value=0).to_markdown())
    if "mp_verified" in e: out.append("\nball-arithmetic sample (mp_verified): %s" % e[e.mp_verified.notna()].mp_verified.value_counts().to_dict())
if len(f) and len(a):
    sec("Arm (f) exploratory: exact q(1..3) as equalities, 32/node, global test")
    af = a[a.wpn == 32]; rows = []
    for cell, g in f.groupby("cell"):
        g0 = af[af.cell == cell]; rows.append(dict(cell=cell, fam=g.fam.iloc[0], n=len(g), DROP_with_q123=int((g.g_outcome == "DROP").sum()), DROP_without=int((g0.g_outcome == "DROP").sum())))
    out.append(pd.DataFrame(rows).to_markdown(index=False))
open(os.path.join(RES, "BD_tables.md"), "w").write("\n".join(out)); print("\n".join(out))
