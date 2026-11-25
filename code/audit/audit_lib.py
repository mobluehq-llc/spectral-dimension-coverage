"""P22 real-data audit: generic adapter that applies E0 / E1q / E4c / E6 (unchanged p22 code) to ANY log-curve with a bootstrap.

Convention: the p22 estimators estimate d = -2 * dy/dx for y = log P.  Here we pass y = -0.5*log F (F = MSD or DFA fluctuation)
so that d = d log F / d log s = alpha (MSD exponent) or H (DFA slope).  x = log scale.  No estimator code is modified; only the
equivalence margin DELTA is set per call (p22 modules import the constant by value, so it is patched in each module).
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); P22 = os.path.dirname(HERE)
sys.path.insert(0, P22)
import numpy as np
import p22lib as L
import p22amend as A
import p22shape as S
import p22honest as H
import p22amend5 as A5
import p22amend6 as A6

H_QUAD = A.H_QUAD


def set_delta(v):
    for m in (L, S, H, A5, A6):
        m.DELTA = v


def wts(Ys):
    v = np.maximum(Ys.var(0, ddof=1), 1e-12)
    w = 1.0 / v
    return w / w.mean()


def _bands_from_score(r, dhat, Ds):
    se = np.maximum(Ds.std(0, ddof=1), 1e-9)
    c = r_crit = np.quantile((np.abs(Ds - dhat) / se).max(1), 0.95)
    return dhat - c * se, dhat + c * se


def run_E1q(x, yA, YsA, yB, YsB, truth=None):
    wA = wts(YsA); best = None
    for h in H_QUAD:
        val, _ = A.poly_rows(x, wA, h, 2); hd = np.diag(val)
        if np.any(np.abs(1 - hd) < 1e-8): continue
        cv = (wA * ((yA - val @ yA) / (1 - hd)) ** 2).sum()
        if best is None or cv < best[0]: best = (cv, h)
    if best is None: return L.abst("E1q", "cvfail"), None
    h = best[1]; _, sl = A.poly_rows(x, wA, h, 2)
    dhat = -2 * sl @ yB; Ds = -2 * YsB @ sl.T
    r = L.score(dhat, Ds, truth, "E1q"); r["h"] = h
    lo, hi = _bands_from_score(r, dhat, Ds)
    return r, dict(d=dhat, lo=lo, hi=hi, x=x)


def run_E0(x, yAll, YsAll, truth_segments=None):
    J = len(x); w = wts(YsAll)
    def lsrow(idx):
        xi, wi = x[idx], w[idx]; xb = (wi * xi).sum() / wi.sum()
        row = np.zeros(J); row[idx] = wi * (xi - xb) / (wi * (xi - xb) ** 2).sum()
        return row, xi.mean()
    half = (J + 1) // 2
    r1, x1 = lsrow(np.arange(half)); r2, x2 = lsrow(np.arange(half, J))
    R = np.stack([r1, r2]); dh = -2 * R @ yAll; Ds = -2 * YsAll @ R.T
    r = L.score(dh, Ds, truth_segments, "E0")
    r["x_seg"] = ";".join("%.3f" % v for v in (x1, x2))
    se = np.maximum(Ds.std(0, ddof=1), 1e-9); c = np.quantile((np.abs(Ds - dh) / se).max(1), 0.95)
    return r, dict(d=dh, lo=dh - c * se, hi=dh + c * se, x=np.array([x1, x2]))


def run_E4c(x, yA, YsA, yB, YsB, truth_d=None, B_unused=None):
    J = len(x); wA = wts(YsA)
    M = H.M_heldout(x, yA, wA)
    if M is None: return L.abst("E4c", "cvfail"), None
    Sm = np.cov(YsB, rowvar=False)
    W = np.zeros((J, J)); sd = np.zeros(J); Bj = np.zeros(J)
    for j in range(J):
        rr = H.honest_weights(x, x[j], Sm, M)
        if rr is None: return L.abst("E4c", "solver"), None
        W[j], sd[j], Bj[j] = rr
    sd = np.maximum(sd, 1e-12); bh = W @ yB; Bs = YsB @ W.T
    c = np.quantile((np.abs(Bs - bh) / sd).max(1), 0.95); hwj = c * sd + Bj
    r = dict(est="E4c", abstain=0, J=J, hw=float(2 * hwj.mean()), se_mean=float(2 * sd.mean()), bias_bound=float(2 * Bj.mean()), M=float(M), crit=float(c))
    Cc = W[1:] - W[[0]]; sdc = np.maximum(np.sqrt(np.einsum("ij,jk,ik->i", Cc, Sm, Cc)), 1e-12); Bc = Bj[1:] + Bj[0]
    dc = Cc @ yB; Dc = YsB @ Cc.T; cc = np.quantile((np.abs(Dc - dc) / sdc).max(1), 0.95); hc = cc * sdc + Bc
    dl = -2 * dc; hcd = 2 * hc; lo, hi = dl - hcd, dl + hcd
    r.update(decl_neg=int((hi < -L.DELTA).any()), decl_pos=int((lo > L.DELTA).any()), decl_zero=int(((hi < 0) | (lo > 0)).any()),
             dmin=float(dl.min()), dmax=float(dl.max()), dlo_min=float(lo.min()), dhi_max=float(hi.max()), jmin_lo=float(lo[np.argmin(dl)]), jmax_hi=float(hi[np.argmax(dl)]))
    d = -2 * bh; dlo = d - 2 * hwj; dhi = d + 2 * hwj
    if truth_d is not None:
        r["cov_simult"] = int(((truth_d >= dlo) & (truth_d <= dhi)).all()); r["cov_point"] = float(((truth_d >= dlo) & (truth_d <= dhi)).mean())
        dt = truth_d[1:] - truth_d[0]; r["delta_cov"] = int(((dt >= lo) & (dt <= hi)).all())
    return r, dict(d=d, lo=dlo, hi=dhi, x=x)


def run_E6(x, yA, YsA, yB, YsB, truth_d=None):
    SA = np.cov(YsA, rowvar=False)
    r, cLo, cHi, SAi = A6.e6_bands(x, yA, SA, yB, YsB, truth_d)
    r.update(e6_location(x, yA, SAi, yB, YsB))
    return r, dict(d=0.5 * (cLo + cHi), lo=cLo, hi=cHi, x=x)


def e6_location(x, yA, SAi, yB, YsB):
    """Logistic-step location m (log scale) over the confidence set (RSS within 99% chi2_2 of the best fit on A), with step size c and the best fit."""
    J = len(x); base = np.stack([np.ones(J), x], 1)
    ms = np.linspace(x[0] - 1.0, x[-1] + 1.0, S.M_GRID_N); fits = []
    SB = np.cov(YsB, rowvar=False); SBi = np.linalg.inv(SB + 1e-10 * np.trace(SB) / J * np.eye(J))
    for m in ms:
        for w in S.W_GRID:
            phi = 0.5 * w * S.softplus((x - m) / w); X = np.column_stack([base, phi])
            K = np.linalg.solve(X.T @ SAi @ X, X.T @ SAi); res = yA - X @ (K @ yA)
            fits.append((m, w, float(res @ SAi @ res)))
    rss = np.array([f[2] for f in fits]); ib = int(np.argmin(rss)); keep = rss - rss.min() <= S.CHI2_2_99
    m0, w0, _ = fits[ib]
    phi = 0.5 * w0 * S.softplus((x - m0) / w0); X = np.column_stack([base, phi])
    th = np.linalg.solve(X.T @ SBi @ X, X.T @ SBi @ yB)
    mk = np.array([f[0] for f, k in zip(fits, keep) if k])
    return dict(e6_m_best=float(m0), e6_w_best=float(w0), e6_step_dB=float(-th[2]), e6_m_lo=float(mk.min()), e6_m_hi=float(mk.max()),
                e6_d_first=float(-2 * th[1] - 0 * th[2]), e6_d_last=float(-2 * th[1] - th[2]))


def run_all(x, A_, B_, All_, truth_d=None, delta=0.25, tags=("E0", "E1q", "E4c", "E6")):
    """A_, B_, All_: (y, Ys) tuples. truth_d: array len J (true local exponent at each scale) or None.
    Returns dict tag -> (row, curve)."""
    set_delta(delta); out = {}
    (yA, YsA), (yB, YsB), (yAll, YsAll) = A_, B_, All_
    if "E0" in tags:
        ts = None
        if truth_d is not None:
            xs = None
        out["E0"] = run_E0(x, yAll, YsAll, None)
        if truth_d is not None:
            r, c = out["E0"]; xseg = c["x"]; td = np.interp(xseg, x, truth_d)
            out["E0"] = run_E0(x, yAll, YsAll, td)
    if "E1q" in tags: out["E1q"] = run_E1q(x, yA, YsA, yB, YsB, truth_d)
    if "E4c" in tags: out["E4c"] = run_E4c(x, yA, YsA, yB, YsB, truth_d)
    if "E6" in tags: out["E6"] = run_E6(x, yA, YsA, yB, YsB, truth_d)
    return out


def summarize(row, with_cov=True):
    keys = ["est", "abstain", "J", "hw", "decl_neg", "decl_pos", "decl_zero", "dmin", "dmax", "cov_simult", "cov_point", "delta_cov", "lof_alarm", "M", "h",
            "e6_m_best", "e6_m_lo", "e6_m_hi", "e6_step_dB", "x_seg"]
    return {k: row.get(k, np.nan) for k in keys}


# ------------------------------------------------------------ minimum-size recipe (Amendment 6), frozen constants read from results
def recipe_params():
    import pandas as pd
    p = pd.read_csv(os.path.join(P22, "results", "amend6_recipe_params.csv"))
    return {r.est: r for _, r in p.iterrows() if r.status == "ok"}


def recipe_view(est, J, S, delta=0.25):
    """Predicted: hw_law(J), ratio r=(S-delta)/hw_law, detect iff J>=Jmin and r>=T.  Also J needed for this S (None if never)."""
    P = recipe_params()
    if est not in P: return dict(hw_law=np.nan, ratio=np.nan, detect=None, J_needed=None, S_min_at_J=np.nan)
    p = P[est]
    hw = lambda j: p.a * (j / 7.0) ** (-p.delta)
    r = (S - delta) / hw(J)
    Jn = None
    for j in range(int(p.Jmin), 400):
        if (S - delta) / hw(j) >= p["T"]: Jn = j; break
    return dict(hw_law=float(hw(J)), ratio=float(r), detect=bool(J >= p.Jmin and r >= p["T"]), J_needed=Jn, S_min_at_J=float(delta + p["T"] * hw(J)), T=float(p["T"]), Jmin=int(p.Jmin))
