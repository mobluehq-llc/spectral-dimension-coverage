#!/usr/bin/env python3
"""Round-3 response, item 5: native-target numbers for E0 and E2 (original replication data, results/raw.csv, wpn 32).
No estimator is re-run.  Usage: round3_native_targets.py --out DIR

Native targets (stated exactly, as the estimators compute them):
  E2  secant on the grid pair (t_k, t_{k+2}):   d_hat_k = -2 ln[q_hat(t_{k+2})/q_hat(t_k)] / ln(t_{k+2}/t_k).
      Its exact counterpart  D2_k = -2 ln[q(t_{k+2})/q(t_k)] / ln(t_{k+2}/t_k)  is, because d(t) = -2 dln q/du (u = ln t), EXACTLY the
      average of d over the log-time interval [ln t_k, ln t_{k+2}] with respect to du.  The manuscript scores E2 against d at the geometric midpoint.
  E0  two weighted least-squares slopes of ln q_hat on u = ln t over the lower and upper half of the window (weights = return counts).
      Exact counterpart D0_s = -2 * WLS slope of ln q(t) on u over the same half, with the counts replaced by their expectation
      w_k = 32 N q(t_k).  (The stored 'cov_pop' columns use the dataset's own observed counts as weights; the two differ negligibly.)
      The manuscript scores E0 against d at exp(mean u) of each half.
Stored per-dataset columns used: cov_simult / cov_point = pointwise target; cov_pop_simult / cov_pop_point = native target;
decl_neg = a drop declared by the estimator's own change test (margin 0.25); abstain = frozen J<6 rule."""
import os, sys
import numpy as np, pandas as pd
from scipy.stats import beta
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
args = sys.argv[1:]; OUT = args[args.index("--out") + 1]; os.makedirs(OUT, exist_ok=True)
RES = os.path.join(HERE, "results")


def cp(k, n):
    k, n = int(k), int(n)
    if n == 0: return (np.nan, np.nan)
    return (beta.ppf(.025, k, n - k + 1) if k > 0 else 0.0, beta.ppf(.975, k + 1, n - k) if k < n else 1.0)


def exact_native(c, J, est, wpn=32):
    """native-target curve, pointwise curve at the matching reference time, for window length J on the cell's grid c['ts_w']."""
    ts = c["ts_w"][:J].astype(float); g = c["graph"]; q, d = g.exact(ts); N = c["N"]
    if est == "E2":
        step = 2; Jd = J - step
        lr = np.log(ts[step:Jd + step] / ts[:Jd]); nat = -2 * np.log(q[step:Jd + step] / q[:Jd]) / lr
        tref = np.sqrt(ts[:Jd] * ts[step:Jd + step]); pw = g.exact(tref)[1]
        return nat, pw
    x = np.log(ts); w = wpn * N * q; half = (J + 1) // 2
    def seg(idx):
        xi, wi = x[idx], w[idx]; xb = (wi * xi).sum() / wi.sum(); sl = (wi * (xi - xb) * np.log(q[idx])).sum() / (wi * (xi - xb) ** 2).sum()
        return -2 * sl, g.exact(np.exp([xi.mean()]))[1][0]
    a = seg(np.arange(half)); b = seg(np.arange(half, J))
    return np.array([a[0], b[0]]), np.array([a[1], b[1]])


def main():
    import run_p22 as R
    raw = pd.read_csv(os.path.join(RES, "raw.csv"), low_memory=False); meta = pd.read_csv(os.path.join(RES, "cells_meta.csv")).set_index("id")
    raw = raw[raw.est.isin(["E0", "E2"])]
    cells = [int(i) for i in sorted(raw.cell.unique())]
    rr = np.random.default_rng(2200); specs = {s["id"]: s for s in R.cell_specs()}
    rows = []
    for cid in cells:
        m = meta.loc[cid]
        if m.kind not in ("crossover", "control") or (isinstance(m.note, str) and m.note != ""): continue
        c = specs[cid]; R.prep(c, rr)
        for est in ("E0", "E2"):
            d = raw[(raw.cell == cid) & (raw.est == est)]
            if len(d) == 0: continue
            n = len(d); ab = int((d.abstain == 1).sum()); a = d[d.abstain == 0]; na = len(a)
            kc = int((a.cov_simult == 1).sum()); kn = int((a.cov_pop_simult == 1).sum()); kd = int((d.decl_neg.fillna(0) == 1).sum()); kp = int((d.decl_pos.fillna(0) == 1).sum())
            J = int(a.J.median()) if na else 0
            if est == "E2": Jwin = J + 2
            else:   # E0 stores J=2 (number of segments); the window length is the expected-count window the estimator would use
                import p22lib as L
                qg = c["graph"].exact(c["ts_w"].astype(float))[0]; Jwin = int(L.window_len(32 * c["N"] * qg, qg, c["N"], c["thr"])); J = Jwin
            r = dict(cell=cid, fam=m.fam, L=int(m.L), kind=m.kind, real=bool(m.real), est=est, n=n, abstain=ab, answered=na, J_median=J,
                     k_cov_pointwise=kc, k_cov_native=kn, k_drop=kd, k_rise=kp)
            for nm, k, den in (("cov_pointwise", kc, na), ("cov_native", kn, na), ("power_all", kd, n), ("power_ans", kd, na)):
                lo, hi = cp(k, den); r[nm] = k / den if den else np.nan; r[nm + "_lo"] = lo; r[nm + "_hi"] = hi
            r["max_drop_pointwise_plan"] = float(m.max_drop)
            if na and Jwin >= 6:
                nat, pw = exact_native(c, Jwin, est)
                r["native_first"] = float(nat[0]); r["native_last"] = float(nat[-1]); r["native_min_change"] = float((nat[1:] - nat[0]).min())
                r["pointwise_first"] = float(pw[0]); r["pointwise_last"] = float(pw[-1]); r["pointwise_min_change"] = float((pw[1:] - pw[0]).min())
                r["native_change_exceeds_margin"] = bool(r["native_min_change"] < -0.25)
            rows.append(r)
        print(cid, m.fam, m.L, flush=True)
    T = pd.DataFrame(rows); T.to_csv(os.path.join(OUT, "r3_native_E0_E2_all_cells.csv"), index=False)
    # bracket the reported thresholds: E0 (2D 192, 3D 96); E2 (2D 192, 3D 128): all tested L within one grid step either side
    thr = {("E0", "aniso2_w16"): 192, ("E0", "aniso3_w8"): 96, ("E2", "aniso2_w16"): 192, ("E2", "aniso3_w8"): 128}
    sel = []
    for (est, fam), L0 in thr.items():
        t = T[(T.est == est) & (T.fam == fam) & (T.kind == "crossover")].sort_values("L")
        Ls = list(t.L); i = Ls.index(L0); keep = Ls[max(i - 2, 0): i + 3]
        x = t[t.L.isin(keep)].copy(); x["reported_threshold_L"] = L0; sel.append(x)
    pd.concat(sel).to_csv(os.path.join(OUT, "r3_native_E0_E2_threshold_brackets.csv"), index=False)


if __name__ == "__main__": main()
