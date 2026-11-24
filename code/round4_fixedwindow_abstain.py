#!/usr/bin/env python3
"""ROUND 4 EXPLORATORY check 3 (spec section 3): fixed-window E4 abstentions in raw_amend10c.csv and bar 10c-B regenerated three ways.
usage: round4_fixedwindow_abstain.py --raw raw_amend10c.csv --meta cells_meta_amend10c.csv --manifest manifest_amend10c.json --out DIR"""
import os, sys, json, hashlib
import numpy as np, pandas as pd
from scipy.stats import beta
a = sys.argv[1:]; g = lambda k: a[a.index(k) + 1]
raw, meta, man, out = g("--raw"), g("--meta"), g("--manifest"), g("--out"); os.makedirs(out, exist_ok=True)
h = hashlib.sha256(open(raw, "rb").read()).hexdigest(); want = json.load(open(man))["csv_sha256"]; assert h == want, (h, want)
LN = []
def P(s=""): print(s); LN.append(str(s))
def cp(k, n):
    return (beta.ppf(.025, k, n - k + 1) if k > 0 else 0.0, beta.ppf(.975, k + 1, n - k) if k < n else 1.0)
def fm(k, n):
    if n == 0: return "0/0"
    lo, hi = cp(k, n); return "%d/%d = %.1f%% (%.1f-%.1f)" % (k, n, 100 * k / n, 100 * lo, 100 * hi)
P(f"EXPLORATORY. raw_amend10c.csv sha256 {h[:16]} == manifest. Realised cells only (cell >= 9000), mode == fixed.")
d = pd.read_csv(raw, low_memory=False); m = pd.read_csv(meta)[["id", "fam", "L", "exact_J", "J_fixed"]].rename(columns={"id": "cell", "J_fixed": "J_fixed_meta"})
d = d[d.cell >= 9000].merge(m, on="cell", how="left"); f = d[d["mode"] == "fixed"].copy()
f["eligible"] = f.J_fixed_meta >= 6                      # exact fixed window has >= 6 points: the frozen J<6 rule does NOT force an abstention
P(f"fixed-window rows by est: {f.est.value_counts().to_dict()}; columns carrying a reason: why={f.why.notna().sum()} non-null")
for est in ("E1q", "E4c", "E4r"):
    x = f[f.est == est]; ab = x[x.abstain == 1]
    P(f"\n== {est}: rows {len(x)}; abstentions {len(ab)} = structural (J_fixed<6) {int((~ab.eligible).sum())} + eligible (J_fixed>=6) {int(ab.eligible.sum())}")
    P("  recorded 'why' among eligible abstentions: " + str(ab[ab.eligible].why.fillna("<none recorded>").value_counts().to_dict()))
    P("  recorded 'why' among structural abstentions: " + str(ab[~ab.eligible].why.fillna("<none recorded>").value_counts().to_dict()))
    P("  'J' column among eligible abstentions: " + str(ab[ab.eligible].J.fillna(-9).value_counts().to_dict()))
    t = ab[ab.eligible].groupby(["fam", "L", "wpn", "J_fixed_meta"]).size().reset_index(name="n_abstain")
    P("  eligible abstentions by cell:\n" + (t.to_string(index=False) if len(t) else "  none"))
    t.to_csv(os.path.join(out, f"round4_fixed_eligible_abstain_{est}.csv"), index=False)
    ab[ab.eligible][["cell", "fam", "L", "wpn", "ds", "J", "J_fixed", "why", "abstain"]].to_csv(os.path.join(out, f"round4_fixed_eligible_rows_{est}.csv"), index=False)
# extra look at the rows of eligible E4c abstentions: what else does the stored row say? (columns that are non-null)
e = f[(f.est == "E4c") & (f.abstain == 1) & f.eligible]
if len(e):
    nn = {c: int(e[c].notna().sum()) for c in e.columns if e[c].notna().any()}; P("\nE4c eligible-abstention rows: non-null columns: " + str(nn))
    P("  distinct J values / J_fixed in those rows: " + str(e.J.value_counts(dropna=False).to_dict()) + " / " + str(e.J_fixed.value_counts(dropna=False).to_dict()))
# bar 10c-B three ways
def pool(est, w, fam, how):
    x = f[(f.est == est) & (f.wpn == w)]
    if fam: x = x[x.fam == fam]
    cov = x[x.abstain == 0].cov_simult.sum(); ans = int((x.abstain == 0).sum())
    if how == "published": return int(cov), ans
    if how == "eligible_miss": return int(cov), ans + int(((x.abstain == 1) & x.eligible).sum())
    if how == "all_miss": return int(cov), len(x)
P("\n== Bar 10c-B: E4c simultaneous coverage >= 93% at each wpn, pooled and per family (point estimate; CP shown; CP-lower strict flag)")
rows = []
for est in ("E4c", "E4r"):
    for fam in (None, "aniso2_w16", "aniso3_w8"):
        for w in (8, 32, 128):
            r = dict(est=est, fam=fam or "pooled", wpn=w)
            for how in ("published", "eligible_miss", "all_miss"):
                k, n = pool(est, w, fam, how); r[how] = fm(k, n); r[how + "_pass"] = bool(n and k / n >= .93 - 1e-12); r[how + "_cpl"] = bool(n and cp(k, n)[0] >= .93 - 1e-12)
            rows.append(r)
T = pd.DataFrame(rows); T.to_csv(os.path.join(out, "round4_bar10cB_three_ways.csv"), index=False)
for est in ("E4c", "E4r"):
    for how in ("published", "eligible_miss", "all_miss"):
        t = T[T.est == est]; P(f"{est} bar B [{how}]: all 9 rows >= 93% (point): {bool(t[how + '_pass'].all())}; also CP-lower >= 93%: {bool(t[how + '_cpl'].all())}; failing rows: " + str([(r.fam, r.wpn, getattr(r, how)) for r in t.itertuples() if not getattr(r, how + '_pass')]))
P("\nE4c detail, highest budget (wpn 128), pooled and per family:")
for r in T[(T.est == "E4c") & (T.wpn == 128)].itertuples(): P(f"  {r.fam}: published {r.published}; eligible-failures-as-misses {r.eligible_miss}; all-abstentions-as-misses {r.all_miss}")
open(os.path.join(out, "round4_fixedwindow_report.txt"), "w").write("\n".join(LN))
