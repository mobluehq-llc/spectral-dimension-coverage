"""P22 Amendment 12: Gaspar-Kukal (Physica A 602, 127624, 2022) weighted log-log estimator of the spectral dimension, as published (eqs 17, 23-24),
adapted only where the data differ (a geometric time grid and a data-driven window; see PLAN.md Amendment 12 and GASPAR_KUKAL_COMPARISON.md).
Arms (all on the common window rule of p22lib.window_len; abstain if J < 6):
  GKW   : WSSQ, weights w_k = p/(1-p), OLS-type analytic standard error (residual-scaled weighted least squares), band +-1.96 se.   [their estimator + their SE]
  GKL   : LSSQ (their classical reference), unweighted, analytic residual-scaled se.
  GKWm  : WSSQ with the model-based (fixed-scale) se: var(ln p_k) = (1-p_k)/(n p_k), n = walkers (their binomial argument), no residual scaling.
  GKWb  : WSSQ point estimate, block-bootstrap se / max-t band (p22lib.score), weights recomputed on every resample.
  GK2W, GK2L, GK2Wb : two-segment extension (lower / upper half of the window, as E0) to obtain a change test; NOT in their paper.
Truth: exact d at the log-mean time of each segment (as E0g/E0). Also the nominal constant dimension on controls (cov_nom)."""
import numpy as np
import p22lib as L

NOMINAL = {"ring": 1.0, "torus2": 2.0, "torus3": 3.0}


def _fit(x, y, w):
    """weighted LS of y on x over the last axis; returns slope, and (sum w, Sxx_w, weighted residual SS)."""
    sw = w.sum(-1); xb = (w * x).sum(-1) / sw; yb = (w * y).sum(-1) / sw
    dx = x - xb[..., None]
    sxx = (w * dx ** 2).sum(-1)
    slope = (w * dx * (y - yb[..., None])).sum(-1) / sxx
    res = (y - yb[..., None]) - slope[..., None] * dx
    return slope, sxx, (w * res ** 2).sum(-1)


def _weights(p, kind):
    return p / (1.0 - p) if kind == "W" else np.ones_like(p)


def _analytic(idxs, x, y, p, kind, ntot, scaled):
    """returns d_hat, se(d_hat), tref (log-mean time) per segment. se(d) = 2 se(slope)."""
    d, se, xr = [], [], []
    for idx in idxs:
        xi, yi, pi = x[idx], y[idx], p[idx]
        w = _weights(pi, kind)
        s, sxx, rss = _fit(xi, yi, w)
        if scaled:
            var = (rss / (len(idx) - 2)) / sxx            # residual-scaled WLS, df = n - 2
        else:
            var = (1.0 / ntot) / sxx                        # w is proportional to the inverse variance n p/(1-p): var(slope) = 1/(n sum w dx^2)
        d.append(-2 * s); se.append(2 * np.sqrt(var)); xr.append(xi.mean())
    return np.array(d), np.array(se), np.exp(np.array(xr))


def _row(tag, d, se, truth_d, nominal, two):
    err = d - truth_d
    z = 1.96 * se
    r = dict(est=tag, abstain=0, J=0, hw=float(z.mean()), se_mean=float(se.mean()), cov_simult=int((np.abs(err) <= z).all()), cov_point=float((np.abs(err) <= z).mean()),
             bias=float(err.mean()), mae=float(np.abs(err).mean()))
    if nominal is not None:
        r["cov_nom"] = int((np.abs(d - nominal) <= z).all()); r["bias_nom"] = float((d - nominal).mean())
    if two:
        dl = d[1] - d[0]; sc = 1.96 * np.sqrt(se[0] ** 2 + se[1] ** 2); lo, hi = dl - sc, dl + sc
        dt = truth_d[1] - truth_d[0]
        r.update(decl_neg=int(hi < -L.DELTA), decl_pos=int(lo > L.DELTA), decl_zero=int((hi < 0) or (lo > 0)), dmin=float(dl), dmax=float(dl),
                 delta_cov=int(lo <= dt <= hi), true_dmin=float(dt))
    return r


def est_gk(Ublk, mwalk, ts, N, rng, truth, B=200, thr=L.THR, nominal=None):
    p = L.grand_mean(Ublk)
    W = sum(u.shape[0] for u in Ublk) * mwalk
    cnt = p * W
    J = L.window_len(cnt, p, N, thr)
    tags = ("GKW", "GKL", "GKWm", "GKWb", "GK2W", "GK2L", "GK2Wb")
    if J < 6:
        return [L.abst(t, "window") for t in tags]
    x = np.log(ts[:J].astype(float)); floor = 0.5 / W
    pj = np.maximum(p[:J], floor); y = np.log(pj)
    half = (J + 1) // 2
    G = [np.arange(J)]; S2 = [np.arange(half), np.arange(half, J)]
    out = []
    for tag, kind, scaled, segs in (("GKW", "W", True, G), ("GKL", "L", True, G), ("GKWm", "W", False, G), ("GK2W", "W", True, S2), ("GK2L", "L", True, S2)):
        d, se, tref = _analytic(segs, x, y, pj, kind, W, scaled)
        out.append(_row(tag, d, se, truth["d_exact"](tref), nominal, len(segs) == 2))
    # bootstrap SE arms (weights recomputed on each resample)
    Ps = np.maximum(L.boot_means(Ublk, B, rng)[:, :J], floor); Ys = np.log(Ps)
    for tag, segs in (("GKWb", G), ("GK2Wb", S2)):
        dh, Dst, trefs = [], [], []
        for idx in segs:
            xi = x[idx]
            s, _, _ = _fit(xi, y[idx], _weights(pj[idx], "W")); dh.append(-2 * s)
            sb, _, _ = _fit(xi[None, :], Ys[:, idx], _weights(Ps[:, idx], "W")); Dst.append(-2 * sb); trefs.append(np.exp(xi.mean()))
        dh = np.array(dh); Dst = np.stack(Dst, 1); trefs = np.array(trefs)
        r = L.score(dh, Dst, truth["d_exact"](trefs), tag)
        if nominal is not None:
            se_b = np.maximum(Dst.std(0, ddof=1), 1e-9); c = np.quantile((np.abs(Dst - dh) / se_b).max(1), 0.95)
            r["cov_nom"] = int((np.abs(dh - nominal) <= c * se_b).all()); r["bias_nom"] = float((dh - nominal).mean())
        out.append(r)
    for r in out:
        r["J"] = J
    return sorted(out, key=lambda r: tags.index(r["est"]))
