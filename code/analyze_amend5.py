#!/usr/bin/env python3
"""P22 Amendment 5 analysis. Reads raw*.csv; writes results/amend5_*.csv and amend5_tables.md."""
import os
import numpy as np, pandas as pd
from p22lib import cp
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
m_old = pd.read_csv(os.path.join(RES, "cells_meta.csv"), keep_default_na=False); m_old["note"] = m_old["note"].astype(str)
m_old = m_old[(m_old.note == "") & m_old.kind.isin(["control", "crossover"])][["id", "fam", "L", "N", "kind", "max_drop", "max_rise"]]
m_new = pd.read_csv(os.path.join(RES, "cells_meta_amend5.csv"))[["id", "fam", "L", "N", "kind", "max_drop", "max_rise"]]
meta = pd.concat([m_old, m_new[m_new.id >= 100]], ignore_index=True)
for c_ in ("max_drop", "max_rise"): meta[c_] = pd.to_numeric(meta[c_], errors="coerce")
def grp(r):
    if r.kind == "control": return "const"
    if r.kind in ("twostep",): return "twostep"
    if r.kind == "brush": return "brush"
    if r.kind == "gasket": return "gasket"
    return "xover_real" if r.max_drop >= 0.6 else "xover_notreal"
meta["group"] = meta.apply(grp, axis=1)
meta["real_drop"] = meta.max_drop >= 0.6; meta["real_rise"] = meta.max_rise >= 0.6
raw5 = pd.read_csv(os.path.join(RES, "raw_amend5.csv"))
raw0 = pd.read_csv(os.path.join(RES, "raw.csv")); raw0 = raw0[raw0.est == "E0"]
raw1 = pd.read_csv(os.path.join(RES, "raw_amend.csv")); raw1 = raw1[raw1.est == "E1q"]
raw3 = pd.read_csv(os.path.join(RES, "raw_amend3.csv")); raw3 = raw3[raw3.est == "E4c"]
old_ids = set(m_old.id)
raw = pd.concat([raw5, raw0[raw0.cell.isin(old_ids)], raw1[raw1.cell.isin(old_ids)], raw3[raw3.cell.isin(old_ids)]], ignore_index=True)
d = raw.merge(meta, left_on="cell", right_on="id")
d["decl_any"] = ((d.decl_neg.fillna(0) + d.decl_pos.fillna(0)) > 0).astype(int)
d["correct"] = np.where(d.real_rise & ~d.real_drop, d.decl_pos, d.decl_neg).astype(float)
d["wrong"] = np.where(d.real_rise & ~d.real_drop, d.decl_neg, d.decl_pos).astype(float)
EST = ["E0", "E1q", "E4c", "E6", "M", "Mb", "CL"]
def ci(k, n): lo, hi = cp(int(k), int(n)); return "%.1f (%.1f-%.1f)" % (100 * k / n, 100 * lo, 100 * hi) if n else "n/a"
out = []
def md(df):
    df = df.fillna("")
    L_ = ["| " + " | ".join(map(str, df.columns)) + " |", "|" + "---|" * len(df.columns)]
    for _, r in df.iterrows(): L_.append("| " + " | ".join(str(v) for v in r.values) + " |")
    return "\n".join(L_)
def P(s=""): out.append(s)

# ---------- E6 lack of fit
e6 = d[(d.est == "E6")].copy(); e6na = e6[e6.abstain == 0].copy()
e6na["miss"] = 1 - e6na.cov_simult
rows = []
for g, x in e6na.groupby("group"):
    rows.append(dict(group=g, n=len(x), abstain=e6[e6.group == g].abstain.mean(), cov=ci(x.cov_simult.sum(), len(x)), alarm=ci(x.lof_alarm.sum(), len(x)),
                     misses=int(x.miss.sum()), alarm_given_miss=ci(x[x.miss == 1].lof_alarm.sum(), int(x.miss.sum())) if x.miss.sum() else "n/a (0 misses)",
                     silent_miss=int(((x.miss == 1) & (x.lof_alarm == 0)).sum()), loud_miss=int(((x.miss == 1) & (x.lof_alarm == 1)).sum()),
                     lof_ratio_median=float(x.lof_ratio.median())))
T6 = pd.DataFrame(rows); T6.to_csv(os.path.join(RES, "amend5_e6_lof_by_group.csv"), index=False)
P("## E6 coverage and lack-of-fit alarm by cell class (non-abstained datasets)"); P(T6.pipe(md)); P()
# per cell
rows = []
for cid, x in e6na.groupby("cell"):
    mm = meta[meta.id == cid].iloc[0]
    rows.append(dict(cell=f"{mm.fam}_{mm.L}", group=mm.group, n=len(x), abstain=e6[e6.cell == cid].abstain.mean(), cov=x.cov_simult.mean(), cov_point=x.cov_point.mean(), alarm=x.lof_alarm.mean(),
                     lof_ratio_med=x.lof_ratio.median(), decl_neg=x.decl_neg.mean(), decl_pos=x.decl_pos.mean(), hw=x.hw.mean(), max_drop=mm.max_drop, max_rise=mm.max_rise))
C6 = pd.DataFrame(rows); C6.to_csv(os.path.join(RES, "amend5_e6_cells.csv"), index=False)
P("## E6 per cell"); P(C6.round(3).pipe(md)); P()
mis = e6na[e6na.group.isin(["twostep", "brush", "gasket"])]
k = int(mis.miss.sum()); P(f"Pooled misspecified cells (twostep+brush+gasket): n={len(mis)}, misses={k}, alarm|miss = " + (ci(mis[mis.miss == 1].lof_alarm.sum(), k) if k else "n/a") + f"; alarm rate = {ci(mis.lof_alarm.sum(), len(mis))}")
ct = e6na[e6na.group == "const"]; P(f"Constant controls: alarm rate (false alarm) = {ci(ct.lof_alarm.sum(), len(ct))}; misses={int(ct.miss.sum())}"); P()
# guarded E6
rows = []
for g, x in e6na.groupby("group"):
    g_decl = ((x.decl_neg + x.decl_pos) > 0) & (x.lof_alarm == 0)
    mm0 = e6[e6.group == g]
    rows.append(dict(group=g, E6_any_decl=ci(x.decl_any.sum(), len(x)), E6g_any_decl=ci(g_decl.sum(), len(x)), E6_correct=ci(x.correct.sum(), len(x)), E6g_correct=ci(((x.correct == 1) & (x.lof_alarm == 0)).sum(), len(x)), E6_wrong=ci(x.wrong.sum(), len(x))))
P("## E6 vs guarded E6 (declaration suppressed when LoF alarm fires), non-abstained"); P(pd.DataFrame(rows).pipe(md)); P()

# ---------- all-estimator comparison per group
rows = []
for est in EST:
    for g in ("const", "xover_real", "xover_notreal", "twostep", "brush", "gasket"):
        x = d[(d.est == est) & (d.group == g)]
        if not len(x): continue
        ok = x[x.abstain == 0]
        r = dict(est=est, group=g, n_cells=x.cell.nunique(), abstain=round(x.abstain.mean(), 3))
        c = ok[ok.cov_simult.notna()]
        r["cov_simult"] = ci(c.cov_simult.sum(), len(c)) if len(c) else "n/a"
        r["cov_point"] = round(c.cov_point.mean(), 3) if len(c) else np.nan
        if g in ("const", "xover_notreal", "gasket"):
            f = ok[ok.decl_neg.notna()]; r["false_change_any"] = ci(f.decl_any.sum(), len(f)) if len(f) else "n/a"; r["false_neg"] = round(f.decl_neg.mean(), 3) if len(f) else np.nan
        else:
            r["power_correct_sign(abst=miss)"] = round(x.correct.fillna(0).mean(), 3); r["wrong_sign"] = round(ok.wrong.mean(), 3) if len(ok) else np.nan
        r["hw"] = round(ok.hw.mean(), 3) if len(ok) else np.nan
        rows.append(r)
TA = pd.DataFrame(rows); TA.to_csv(os.path.join(RES, "amend5_groups_all_estimators.csv"), index=False)
P("## All estimators by cell class (cov_simult = simultaneous coverage of exact d(t) at each estimator's own reference scales; CP 95%)"); P(TA.pipe(md)); P()

# ---------- pass bar on controls for new arms + min L
rows = []
cell_name = lambda c: f"{meta.set_index('id').loc[c].fam}{meta.set_index('id').loc[c].L}"
for est in EST:
    r = dict(est=est)
    g = d[(d.est == est) & (d.group == "const")]; ok = g[g.abstain == 0]
    r["abstain_ctrl"] = round(g.abstain.mean(), 3)
    c = ok[ok.cov_simult.notna()]
    if len(c):
        k = c.cov_simult.sum(); r["cov_ctrl"] = ci(k, len(c))
        pc = c.groupby("cell").cov_simult.agg(["sum", "count"])
        bad = [cell_name(i) for i, v in pc.iterrows() if v["count"] >= 20 and cp(int(v["sum"]), int(v["count"]))[1] < 0.93]
        r["cov_bad_cells"] = ";".join(bad) or "none"; r["C"] = bool(k / len(c) >= 0.93 and not bad)
        if "cov_simult_plateau" in c and c.cov_simult_plateau.notna().any():
            cp_ = c[c.cov_simult_plateau.notna()]; r["cov_ctrl_vs_plateau"] = ci(cp_.cov_simult_plateau.sum(), len(cp_))
    f = ok[ok.decl_neg.notna()]
    if len(f):
        k = f.decl_any.sum(); r["fpr"] = ci(k, len(f))
        pf = f.groupby("cell").decl_any.agg(["sum", "count"])
        badf = [cell_name(i) for i, v in pf.iterrows() if v["count"] >= 20 and cp(int(v["sum"]), int(v["count"]))[0] > 0.05]
        r["fpr_bad_cells"] = ";".join(badf) or "none"; r["F"] = bool(k / len(f) <= 0.05 and not badf)
    cr = d[(d.est == est) & (d.group == "xover_real") & (d.abstain == 0) & d.cov_simult.notna()]
    r["cov_xover_real"] = ci(cr.cov_simult.sum(), len(cr)) if len(cr) else "n/a"
    pr = d[(d.est == est) & (d.group == "xover_notreal") & (d.abstain == 0) & d.decl_neg.notna()]
    r["premature_any"] = ci(pr.decl_any.sum(), len(pr)) if len(pr) else "n/a"
    for fam in ("aniso2_w16", "aniso3_w8"):
        c2 = d[(d.est == est) & (d.fam == fam) & (d.group == "xover_real")]
        res = "n/a"
        if len(c2) and c2.decl_neg.notna().any():
            t = c2.groupby("L").apply(lambda x: pd.Series(dict(k=x.decl_neg.fillna(0).sum(), n=len(x))), include_groups=False).sort_index()
            pw = (t.k / t.n).values
            def first(cond):
                for i in range(len(t)):
                    if cond[i:].all(): return i
            ip = first(pw >= 0.8)
            res = ("%d" % t.index[ip]) if ip is not None else "not reached (>%d)" % t.index.max()
            r["power_at_Lmax_" + fam] = round(float(pw[-1]), 2)
        r["minL_" + fam] = res
    rows.append(r)
TP = pd.DataFrame(rows); TP["PASS_CF"] = TP.C.fillna(False) & TP.F.fillna(False)
TP.to_csv(os.path.join(RES, "amend5_passbar.csv"), index=False)
P("## Pass bar (controls; C = simult. coverage >= 93% pooled and no cell CP-upper < 93%; F = FPR <= 5%), premature, min L"); P(TP.pipe(md)); P()
# E4c detail
for est in ("M", "Mb", "CL"):
    x = d[(d.est == est) & (d.abstain == 0)]
    P(f"{est}: median window/ks columns n/a here; mean |bias| on const controls = {x[x.group=='const'].bias.abs().mean():.3f}")
open(os.path.join(RES, "amend5_tables.md"), "w").write("\n".join(out))
print("\n".join(out))
