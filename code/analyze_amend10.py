#!/usr/bin/env python3
"""P22 Amendment 10 analysis. Usage: analyze_amend10.py a|b|c.   10a: margin-consistent false-change scoring of existing data. Verdict lines printed literally."""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import beta
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
sys.path.insert(0, HERE)
LINES = []
def P(s=""): print(s); LINES.append(str(s))
def cp(k, n):
    k, n = int(k), int(n)
    if n == 0: return (np.nan, np.nan)
    return (beta.ppf(.025, k, n - k + 1) if k > 0 else 0.0, beta.ppf(.975, k + 1, n - k) if k < n else 1.0)
def fmt(k, n):
    lo, hi = cp(k, n); return "%d/%d = %.1f%% (CP %.1f-%.1f)" % (k, n, 100 * k / max(n, 1), 100 * lo, 100 * hi)

# ============================================================ 10a
NR = [14, 15, 16, 17, 21]                    # original non-realised crossover cells (max_drop < 0.6): aniso2_w16 L=32,48,64,96 ; aniso3_w8 L=32
SRC = {"E3": "raw.csv", "E0": "raw.csv", "E1q": "raw_amend.csv", "CL": "raw_amend5.csv"}
MARG = 0.25

def part_a():
    import run_p22 as R
    import p22lib as L
    specs = {s["id"]: s for s in R.cell_specs() if s["id"] in NR}
    rr = np.random.default_rng(2200)
    cells = {}
    for i, s in specs.items():
        R.prep(s, rr); ts = s["ts_w"]; cells[i] = dict(ts=ts, d=s["d_exact"](ts), dfun=s["d_exact"], fam=s["fam"], L=s["L"], max_drop=s["max_drop"])
    rows = []
    for est, f in SRC.items():
        d = pd.read_csv(os.path.join(RES, f), low_memory=False); d = d[(d.est == est) & d.cell.isin(NR)]
        for _, r in d.iterrows():
            c = cells[int(r.cell)]; ts, dx = c["ts"], c["d"]
            out = dict(est=est, cell=int(r.cell), ds=int(r.ds), abstain=int(r.abstain), fam=c["fam"], L=c["L"], max_drop=c["max_drop"])
            if r.abstain == 1:
                rows.append(out); continue
            dn, dp = int(r.decl_neg), int(r.decl_pos)
            if est in ("E0", "E1q"):
                J = int(r.J); jj = np.arange(J)
                if est == "E0":
                    half = (J + 1) // 2; x = np.log(ts[:J].astype(float))
                    cen = np.exp([x[:half].mean(), x[half:].mean()]); dcen = float(np.diff(c["dfun"](cen))[0])
                else:
                    dcen = float(np.diff(c["dfun"](np.exp(np.log(ts[:J].astype(float))[[0, 1]])))[0]) if False else np.nan
            else:
                tl = r.tsplit ** 2 / 24.0; jj = np.nonzero(ts <= tl * (1 + 1e-9))[0]
                if len(jj) < 2: jj = np.arange(2)
                J = len(jj); tt = ts[jj].astype(float)
                cen = np.exp([np.log(tt[tt <= r.tsplit]).mean(), np.log(tt[tt > r.tsplit]).mean()]); dcen = float(np.diff(c["dfun"](cen))[0])
            exc = dx[jj][1:] - dx[0]                                   # Delta_j = d(t_j)-d(t_0), j=1..J-1, t_j in the estimator's own window
            out.update(J=J, ex_min=float(exc.min()), ex_max=float(exc.max()), ex_cen=dcen, decl_neg=dn, decl_pos=dp)
            rows.append(out)
    D = pd.DataFrame(rows)
    return D, cells

def add_centre(D, cells, ctr):
    return D

def classify(D, mode="window"):
    D = D.copy()
    ex_drop = -D.ex_min; ex_rise = D.ex_max
    D["declared"] = ((D.decl_neg == 1) | (D.decl_pos == 1)).astype(int)
    tn = (D.decl_neg == 1) & (ex_drop > MARG); tp = (D.decl_pos == 1) & (ex_rise > MARG)
    D["true_chg"] = (tn | tp).astype(int)
    tnc = (D.decl_neg == 1) & (-D.ex_cen > MARG); tpc = (D.decl_pos == 1) & (D.ex_cen > MARG)
    D["true_cen"] = (tnc | tpc).astype(int); D["false_cen"] = ((D.declared == 1) & (D.true_cen == 0)).astype(int)
    D["false_chg"] = ((D.declared == 1) & (D.true_chg == 0)).astype(int)
    return D

def rep_a():
    D, cells = part_a()
    D.to_csv(os.path.join(RES, "amend10a_rows.csv"), index=False)
    P("# 10a: margin-consistent false-change scoring (original non-realised crossover cells)")
    P("true change: declared sign matches and exact excursion (d(t_j)-d(t_0)) beyond 0.25 in that direction at some t_j in the estimator's own tested window; false change: declared otherwise (wrong sign counts as false).")
    P("tested window: E1q = ts[:J] (J per dataset); E0 = ts[:J]; E3/CL = grid times <= tlast of the largest graph, tlast = tsplit^2/24 (tsplit = sqrt(24 tlast) by construction). Excursion taken vs the first grid point t_0=24.")
    P("sensitivity: the same scoring with the excursion defined as d(late centre)-d(early centre) (E3/CL: region centres of the collapse window; E0: segment centres).")
    P("exact excursion over the whole collapse window per cell:")
    for i in NR:
        c = cells[i]; ts = c["ts"]; dx = c["d"]
        P("  cell %d %s L=%d  max_drop(Nq>=6 window, plan)=%.4f ; grid %d..%d ; d(t)-d(24) min over grid<=1.05-window: see rows" % (i, c["fam"], c["L"], c["max_drop"], ts[0], ts[-1]))
    DD = classify(D[D.abstain == 0])
    tab = []
    for est in SRC:
        P("\n## %s" % est)
        for i in NR:
            g = D[(D.est == est) & (D.cell == i)]; ga = DD[(DD.est == est) & (DD.cell == i)]
            n, na = len(ga), int(len(g) - len(ga))
            ex = ga.ex_min.median() if len(ga) else np.nan
            P("  cell %2d %-10s L=%-3d drop(plan)=%.3f | abstain %d, non-abstained %d, median exact min-excursion over window %.3f | no-decl %d | TRUE %s | FALSE %s | declared %s" % (
                i, cells[i]["fam"], cells[i]["L"], cells[i]["max_drop"], na, n, ex, int((ga.declared == 0).sum()), fmt(ga.true_chg.sum(), n), fmt(ga.false_chg.sum(), n), fmt(ga.declared.sum(), n)))
            tab.append(dict(est=est, cell=i, fam=cells[i]["fam"], L=cells[i]["L"], max_drop=cells[i]["max_drop"], n_total=len(g), n_abstain=na, n=n, n_nodecl=int((ga.declared == 0).sum()),
                            n_true=int(ga.true_chg.sum()), n_false=int(ga.false_chg.sum()), n_declared=int(ga.declared.sum()), n_false_cen=int(ga.false_cen.sum()), n_true_cen=int(ga.true_cen.sum()), med_ex_cen=float(ga.ex_cen.median()) if len(ga) else np.nan))
    T = pd.DataFrame(tab); T.to_csv(os.path.join(RES, "amend10a_table.csv"), index=False)
    # pooled sets
    P("\n## pooled sets (false-change rate; bar: >=20%% and CP-lower >10%% with >=30 non-abstained)")
    sets = {"near-flat, plan wording 'drop <= 0.02' read as the two cells with drop 0.005 and 0.021 (ids 14,15)": [14, 15],
            "near-flat, strictly drop <= 0.02 (id 14 only)": [14],
            "all five non-realised cells": NR, "drop 0.34 and 0.57 cells (ids 17,21)": [17, 21]}
    verd = {}
    for name, ids in sets.items():
        P("\n" + name)
        for est in SRC:
            x = T[(T.est == est) & T.cell.isin(ids)]; n = int(x.n.sum()); k = int(x.n_false.sum()); kt = int(x.n_true.sum())
            lo, hi = cp(k, n)
            ok = n >= 30 and k / max(n, 1) >= 0.20 and lo > 0.10
            P("  %-4s false %s | true %s | n=%d -> %s" % (est, fmt(k, n), fmt(kt, n), n, "bar met (headline survives for this arm)" if ok else ("bar NOT met" if n >= 30 else "n<30: cannot assess (bar requires >=30)")))
            verd[(name, est)] = (ok, n, k)
            kc = int(x.n_false_cen.sum()); lc, _ = cp(kc, n)
            P("       sensitivity (excursion between the two regions' log-mean centres instead of the window): false %s, true %d" % (fmt(kc, n), int(x.n_true_cen.sum())))
    k1 = list(sets)[0]; k2 = list(sets)[1]
    P("\nVERDICT 10a (plan wording read as ids 14,15): " + "; ".join("%s: %s" % (e, "SURVIVES" if verd[(k1, e)][0] else ("NOT ASSESSABLE n=%d<30 (all abstain)" % verd[(k1, e)][1] if verd[(k1, e)][1] < 30 else "headline NOT SUPPORTED: bar fails (n=%d)" % verd[(k1, e)][1])) for e in SRC))
    P("VERDICT 10a (strict <=0.02, id 14 only): " + "; ".join("%s: %s" % (e, "SURVIVES" if verd[(k2, e)][0] else ("NOT ASSESSABLE n=%d<30 (all abstain)" % verd[(k2, e)][1] if verd[(k2, e)][1] < 30 else "headline NOT SUPPORTED: bar fails (n=%d)" % verd[(k2, e)][1])) for e in SRC))
    open(os.path.join(RES, "amend10a_report.txt"), "w").write("\n".join(LINES))

if __name__ == "__main__":
    which = sys.argv[1]
    if which == "a": rep_a()
