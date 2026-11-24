#!/usr/bin/env python3
"""P22 Amendment 4 analysis: union table of all estimators. Reads raw.csv + raw_amend.csv, writes results/amend_*.csv."""
import os
import numpy as np, pandas as pd
from p22lib import cp
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
meta = pd.read_csv(os.path.join(RES, "cells_meta.csv"), keep_default_na=False); meta["note"] = meta["note"].astype(str)
raw = pd.concat([pd.read_csv(os.path.join(RES, "raw.csv")), pd.read_csv(os.path.join(RES, "raw_amend.csv")), pd.read_csv(os.path.join(RES, "raw_amend3.csv")), pd.read_csv(os.path.join(RES, "raw_amend4.csv"))], ignore_index=True)
raw = raw[raw.est != "E1_dup"]
d = raw.merge(meta[["id", "fam", "L", "N", "kind", "note", "real", "max_drop"]], left_on="cell", right_on="id")
d = d[(d.note == "") & d.kind.isin(["control", "crossover"])].copy()
d["decl"] = ((d.decl_neg.fillna(0) + d.decl_pos.fillna(0)) > 0).astype(int)
EST = ["E0", "E1", "E1x", "E2", "E2h", "E3", "E3s", "E1u", "E1q", "E1qu", "E4o", "E4c", "E4r", "E5qu", "E5q", "E6"]
ctrl = d[d.kind == "control"]; cross = d[d.kind == "crossover"]
def ci(k, n): lo, hi = cp(int(k), int(n)); return "%.1f (%.1f-%.1f)" % (100 * k / n, 100 * lo, 100 * hi) if n else "n/a"
rows, cellrows = [], []
for est in EST:
    r = dict(est=est)
    g = ctrl[(ctrl.est == est)]; ok = g[g.abstain == 0]
    r["abstain_ctrl"] = g.abstain.mean()
    gc = ok[ok.cov_simult.notna()]
    if len(gc):
        k = gc.cov_simult.sum(); r["cov_ctrl"] = ci(k, len(gc)); r["cov_ctrl_pt"] = k / len(gc)
        pc = gc.groupby("cell").cov_simult.agg(["sum", "count"])
        bad = [f"{meta.set_index('id').loc[c].fam}{meta.set_index('id').loc[c].L}" for c, x in pc.iterrows() if x["count"] >= 20 and cp(int(x["sum"]), int(x["count"]))[1] < 0.93]
        r["cov_bad_cells"] = ";".join(bad) or "none"; r["C"] = bool(k / len(gc) >= 0.93 and not bad)
    gf = ok[ok.decl_neg.notna()]
    if len(gf):
        k = gf.decl.sum(); r["fpr"] = ci(k, len(gf)); r["fpr_pt"] = k / len(gf)
        pf = gf.groupby("cell").decl.agg(["sum", "count"])
        badf = [c for c, x in pf.iterrows() if x["count"] >= 20 and cp(int(x["sum"]), int(x["count"]))[0] > 0.05]
        r["fpr_bad_cells"] = ";".join(f"{meta.set_index('id').loc[c].fam}{meta.set_index('id').loc[c].L}" for c in badf) or "none"
        r["F"] = bool(k / len(gf) <= 0.05 and not badf)
    # coverage vs exact d on REALISED crossovers
    cr = cross[(cross.est == est) & (cross.real == True) & (cross.abstain == 0) & cross.cov_simult.notna()]
    if len(cr):
        k = cr.cov_simult.sum(); r["cov_cross_real"] = ci(k, len(cr)); r["cov_cross_real_pt"] = k / len(cr)
        pc = cr.groupby("cell").cov_simult.mean()
        r["cov_cross_range"] = "%.0f-%.0f%%" % (100 * pc.min(), 100 * pc.max())
    else:
        r["cov_cross_real"] = "n/a (target = thermodynamic d)"
    # premature: non-real crossover cells, non-abstained
    pr = cross[(cross.est == est) & (cross.real == False) & (cross.abstain == 0) & cross.decl_neg.notna()]
    if len(pr):
        r["premature_neg"] = ci(pr.decl_neg.sum(), len(pr)); r["premature_neg_pt"] = pr.decl_neg.mean(); r["premature_any_pt"] = pr.decl.mean(); r["n_premature"] = len(pr)
    # power (abstention = miss)
    for fam in ("aniso2_w16", "aniso3_w8"):
        c = cross[(cross.est == est) & (cross.fam == fam) & (cross.real == True)]
        out = "n/a"
        if len(c) and c.decl_neg.notna().any():
            t = c.groupby("L").apply(lambda x: pd.Series(dict(k=x.decl_neg.fillna(0).sum(), n=len(x), N=x.N.iloc[0])), include_groups=False).sort_index()
            pw = (t.k / t.n).values; lo = np.array([cp(int(a), int(b))[0] for a, b in zip(t.k, t.n)])
            def first(cond):
                for i in range(len(t)):
                    if cond[i:].all(): return i
            ip = first(pw >= 0.8); isx = first(lo >= 0.8)
            out = ("%d" % t.index[ip]) if ip is not None else "not reached (>%d)" % t.index.max()
            r["minL_strict_" + fam] = ("%d" % t.index[isx]) if isx is not None else "not reached"
        r["minL_" + fam] = out
    rows.append(r)
T = pd.DataFrame(rows); T["PASS_CF"] = T.C.fillna(False) & T.F.fillna(False)
T.to_csv(os.path.join(RES, "amend4_final_table.csv"), index=False)
# per-cell for the new + hero estimators
cc = []
for (cid, est), g in d.groupby(["cell", "est"]):
    m = meta[meta.id == cid].iloc[0]; ok = g[g.abstain == 0]
    cc.append(dict(fam=m.fam, L=m.L, kind=m.kind, real=m.real, max_drop=m.max_drop, est=est, abstain=g.abstain.mean(),
                   cov=ok.cov_simult.mean() if ok.cov_simult.notna().any() else np.nan,
                   cov_pop=ok.cov_pop_simult.mean() if "cov_pop_simult" in ok and ok.cov_pop_simult.notna().any() else np.nan,
                   decl_neg=ok.decl_neg.mean() if ok.decl_neg.notna().any() else np.nan, decl_any=ok.decl.mean() if len(ok) and ok.decl_neg.notna().any() else np.nan,
                   power=g.decl_neg.fillna(0).mean(), bias=ok.bias.mean() if "bias" in ok else np.nan, hw=ok.hw.mean() if len(ok) else np.nan, h=ok.h.mean() if "h" in ok and ok.h.notna().any() else np.nan))
CC = pd.DataFrame(cc).sort_values(["kind", "fam", "L", "est"]); CC.to_csv(os.path.join(RES, "amend4_cell_summary.csv"), index=False)
pd.set_option("display.width", 300, "display.max_columns", 60, "display.max_colwidth", 60)
print(T.drop(columns=[c for c in T if c.endswith("_pt")]).T.to_string())
print(CC[(CC.est.isin(["E1q", "E1u", "E4o", "E4c", "E4r", "E5qu", "E5q", "E6"])) & (CC.kind == "crossover")].round(3).to_string())
