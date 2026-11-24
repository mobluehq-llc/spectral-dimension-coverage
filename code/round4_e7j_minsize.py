#!/usr/bin/env python3
"""ROUND 4 EXPLORATORY check 4 (spec section 4): what the stored results say about the E7j (3D, aniso3_w8) minimum size. Reads stored files only.
usage: round4_e7j_minsize.py --res results/ --a10 results_amend10_11/ --out DIR"""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import beta
a = sys.argv[1:]; g = lambda k: a[a.index(k) + 1]; RES, A10, OUT = g("--res"), g("--a10"), g("--out"); os.makedirs(OUT, exist_ok=True)
LN = []
def P(s=""): print(s); LN.append(str(s))
def cp(k, n): return (beta.ppf(.025, k, n - k + 1) if k > 0 else 0.0, beta.ppf(.975, k + 1, n - k) if k < n else 1.0)
def fm(k, n):
    if n == 0: return "0/0"
    lo, hi = cp(k, n); return "%d/%d = %.1f%% (%.1f-%.1f)" % (k, n, 100 * k / n, 100 * lo, 100 * hi)
d = pd.read_csv(os.path.join(A10, "raw_amend10b.csv"), low_memory=False)
for c in ("Bj_neg", "Bj_pos", "Bj_amb"): d[c] = d[c].fillna(0)
d["B_ok"] = d.B_status == 0; d["j_drop"] = (d.B_ok & (d.Bj_neg == 1)).astype(int)
x = d[d.fam == "aniso3_w8"]
P("EXPLORATORY. Stored 10b outer-relaxation-B certified E7j drop (the 10b-M scoring: no band = no detection), 3D aniso3_w8, by L and selection group.")
P("10b groups present: " + str(d.group.value_counts().to_dict()))
m9 = pd.read_csv(os.path.join(RES, "cells_meta_amend9_a.csv")); r9 = pd.read_csv(os.path.join(RES, "raw_amend9_a.csv"), low_memory=False)
r9 = r9.merge(m9[["id", "fam", "L"]].rename(columns={"id": "cell"}), on="cell"); r9 = r9[(r9.fam == "aniso3_w8")]
r9["disc_drop"] = (r9.decl_neg == 1).astype(int)
rows = []
for L, t in x.groupby("L"):
    k, n = int(t.j_drop.sum()), len(t)
    r9L = r9[r9.L == L]; full9 = fm(int(r9L.disc_drop.sum()), len(r9L)) if len(r9L) else "n/a"
    # same datasets under the 9a discretised scoring
    mm = t.merge(r9L[["ds", "disc_drop"]], on="ds", how="left"); nm = int(mm.disc_drop.notna().sum()); kd = int(mm.disc_drop.sum())
    rows.append(dict(L=int(L), n_10b=n, groups=str(t.group.value_counts().to_dict()), outerB_drop=fm(k, n), outerB_k=k, same_datasets_discretised_drop=fm(kd, nm) if nm else "n/a", all200_discretised_drop_9a=full9, is_full_200=n == 200))
T = pd.DataFrame(rows); T.to_csv(os.path.join(OUT, "round4_e7j_minsize_table.csv"), index=False)
P("\nL | n in 10b | selection groups | 10b outer-B certified E7j drop (fixed denominator) | same datasets, 9a discretised E7j | all 200 datasets of that L, 9a discretised E7j")
for r in rows: P(f"L={r['L']:>4} | n={r['n_10b']:>3}{'  (FULL 200)' if r['is_full_200'] else '  (SELECTED SUBSET)'} | {r['groups']} | {r['outerB_drop']} | {r['same_datasets_discretised_drop']} | {r['all200_discretised_drop_9a']}")
P("\nBar: certified power >= 0.8 at L=160,192,256 (pooled point rule used by 10b-M).")
for L in (160, 192, 256):
    t = x[x.L == L]; P(f"  L={L}: {fm(int(t.j_drop.sum()), len(t))}  (CP-lower >= 0.8: {cp(int(t.j_drop.sum()), len(t))[0] >= .8})")
P("\nSmaller sizes: the 10b outer-B datasets are NOT a random sample of the 200 (they were chosen as atom-sensitive or withheld ones, plus a seeded stable draw); their rates are not power estimates.")
P("Only the 9a discretised scoring covers all 200 datasets at every size; there every size < 160 is below 0.8 (see last column).")
open(os.path.join(OUT, "round4_e7j_minsize_report.txt"), "w").write("\n".join(LN))
