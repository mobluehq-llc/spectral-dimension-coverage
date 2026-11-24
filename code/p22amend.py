"""P22 Amendment 1 estimators (see PLAN.md AMENDMENT 1): E1u, E1q, E1qu (bias-corrected local polynomial) and E3s (single-size E3)."""
import numpy as np
import p22lib as L
from p22lib import grand_mean, boot_means, score, abst, window_len

H_LIN = L.H_GRID
H_QUAD = (0.7, 1.0, 1.4, 2.0, 2.8, 4.0)


def poly_rows(x, w, h, deg):
    """Weighted local polynomial (Gaussian kernel, weights w) at each x_j. Returns (val rows, slope rows), each J x J, linear in y."""
    J = len(x)
    val = np.zeros((J, J)); sl = np.zeros((J, J))
    for j in range(J):
        z = x - x[j]
        u = np.exp(-0.5 * (z / h) ** 2) * w
        X = np.vander(z, deg + 1, increasing=True)
        A = X.T @ (u[:, None] * X)
        try:
            M = np.linalg.solve(A, X.T * u[None, :])
        except np.linalg.LinAlgError:
            M = np.linalg.pinv(A) @ (X.T * u[None, :])
        val[j] = M[0]; sl[j] = M[1]
    return val, sl


def est_E1b(Ublk, mwalk, ts, N, rng, truth, deg, under, B=200, thr=L.THR, tag="E1q"):
    UA = [u[0::2] for u in Ublk]; UB = [u[1::2] for u in Ublk]
    pA = grand_mean(UA)
    WA = sum(u.shape[0] for u in UA) * mwalk
    cntA = pA * WA
    J = window_len(cntA, pA, N, thr)
    if J < 6:
        return abst(tag, "window")
    x = np.log(ts[:J].astype(float)); yA = np.log(pA[:J]); wA = cntA[:J]
    best = None
    for h in (H_LIN if deg == 1 else H_QUAD):
        val, _ = poly_rows(x, wA, h, deg)
        hd = np.diag(val)
        if np.any(np.abs(1 - hd) < 1e-8):
            continue
        cv = (wA * ((yA - val @ yA) / (1 - hd)) ** 2).sum()
        if best is None or cv < best[0]:
            best = (cv, h)
    if best is None:
        return abst(tag, "cvfail")
    h = best[1] * under
    _, sl = poly_rows(x, wA, h, deg)
    WB = sum(u.shape[0] for u in UB) * mwalk
    floor = 0.5 / WB
    pB = np.maximum(grand_mean(UB), floor)[:J]
    dhat = -2 * sl @ np.log(pB)
    Ps = np.maximum(boot_means(UB, B, rng)[:, :J], floor)
    Ds = -2 * np.log(Ps) @ sl.T
    dtrue = truth["d_exact"](ts[:J])
    dpop = -2 * sl @ np.log(truth["q_exact"](ts[:J]))
    r = score(dhat, Ds, dtrue, tag, dpop)
    r.update(h=h, tlast=int(ts[J - 1]))
    return r


def est_E3s(Ublk, mwalk, ts, N, rng, truth, B=80, thr=1.05, tag="E3s"):
    p = grand_mean(Ublk)
    W = sum(u.shape[0] for u in Ublk) * mwalk
    cnt = p * W
    J = window_len(cnt, p, N, thr)
    if J < 6:
        return abst(tag, "window")
    tsplit = np.sqrt(24.0 * ts[J - 1])
    t = ts[:J].astype(float)
    e = np.nonzero(t <= tsplit)[0]; l = np.nonzero(t > tsplit)[0]
    if len(e) < 3 or len(l) < 3:
        return abst(tag, "halves")
    x = np.log(t); w = cnt[:J]
    def lsrow(idx):
        xi, wi = x[idx], w[idx]
        xb = (wi * xi).sum() / wi.sum()
        row = np.zeros(J); row[idx] = wi * (xi - xb) / (wi * (xi - xb) ** 2).sum()
        return row, xi.mean()
    r1, x1 = lsrow(e); r2, x2 = lsrow(l)
    R = np.stack([r1, r2])
    floor = 0.5 / W
    y = np.log(np.maximum(p[:J], floor))
    Ys = np.log(np.maximum(boot_means(Ublk, B, rng)[:, :J], floor))
    qe = np.log(truth["q_exact"](ts[:J]))
    return score(-2 * (R @ y), -2 * Ys @ R.T, truth["d_exact"](np.exp([x1, x2])), tag, -2 * (R @ qe))
