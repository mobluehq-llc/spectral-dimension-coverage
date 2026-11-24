"""P22 Amendment 4: E5 (monotone multiscale band, Duembgen-style) and E6 (split-sample logistic-step inference). See PLAN.md AMENDMENT 4."""
import numpy as np
import p22lib as L
import p22amend as A
from p22lib import grand_mean, boot_means, abst, window_len, DELTA

M_GRID_N = 21
W_GRID = (0.1, 0.2, 0.4, 0.8, 1.5, 2.5)
CHI2_2_99 = 9.2103
CHI2_2_96 = 6.4378   # chi2_2 quantile 0.96 = -2 ln(0.04)
SHRINK = 0.1


def _intervals(J):
    return [(i, k) for i in range(J) for k in range(i, J)]


def monotone_band(dhat, Ds, decreasing=True, level=0.95):
    """returns lo, hi (per grid point), crit c. Increasing handled by sign flip."""
    if not decreasing:
        lo, hi, c = monotone_band(-dhat, -Ds, True, level)
        return -hi, -lo, c
    J = len(dhat)
    iv = _intervals(J)
    cs = np.concatenate([np.zeros((Ds.shape[0], 1)), np.cumsum(Ds, 1)], 1)
    c0 = np.concatenate([[0.0], np.cumsum(dhat)])
    ia = np.array([a for a, b in iv]); ib = np.array([b for a, b in iv])
    n = (ib - ia + 1).astype(float)
    mB = (cs[:, ib + 1] - cs[:, ia]) / n
    m0 = (c0[ib + 1] - c0[ia]) / n
    sd = np.maximum(mB.std(0, ddof=1), 1e-9)
    c = float(np.quantile((np.abs(mB - m0) / sd).max(1), level))
    up = m0 + c * sd; lw = m0 - c * sd
    U = np.full(J, np.inf); Lo = np.full(J, -np.inf)
    for j in range(J):
        s1 = ib <= j
        if s1.any(): U[j] = up[s1].min()
        s2 = ia >= j
        if s2.any(): Lo[j] = lw[s2].max()
    return Lo, U, c


def e5_core(dhat, Ds, dtrue, tag, decreasing=True):
    Lo, U, c = monotone_band(dhat, Ds, decreasing)
    r = dict(est=tag, abstain=0, J=len(dhat), crit=c, hw=float(np.mean(U - Lo) / 2), se_mean=float(Ds.std(0, ddof=1).mean()))
    if decreasing:
        drop = Lo[0] - U[-1]; r.update(decl_neg=int(drop > DELTA), decl_pos=0, decl_zero=int(drop > 0), dmin=float(-max(drop, 0)), dmax=0.0)
    else:
        rise = Lo[-1] - U[0]; r.update(decl_neg=0, decl_pos=int(rise > DELTA), decl_zero=int(rise > 0), dmin=0.0, dmax=float(max(rise, 0)))
    r["band_drop_lo"] = float(Lo[0] - U[-1]) if decreasing else float(Lo[-1] - U[0])
    if dtrue is not None:
        r["cov_simult"] = int(((dtrue >= Lo) & (dtrue <= U)).all()); r["cov_point"] = float(((dtrue >= Lo) & (dtrue <= U)).mean())
        r["bias"] = float((dhat - dtrue).mean()); r["mae"] = float(np.abs(dhat - dtrue).mean())
        if decreasing:
            dt = dtrue[1:] - dtrue[0]; lo, hi = Lo[1:] - U[0], U[1:] - Lo[0]
        else:
            dt = dtrue[1:] - dtrue[0]; lo, hi = Lo[1:] - U[0], U[1:] - Lo[0]
        r["delta_cov"] = int(((dt >= lo) & (dt <= hi)).all()); r["true_dmin"] = float(dt.min())
    return r


def slope_ops(UA_p, wA, x, thr_ok=True):
    """chosen h_cv (quadratic) on A; returns dict tag->slope operator."""
    best = None
    for h in A.H_QUAD:
        val, _ = A.poly_rows(x, wA, h, 2)
        hd = np.diag(val)
        if np.any(np.abs(1 - hd) < 1e-8): continue
        cv = (wA * ((UA_p - val @ UA_p) / (1 - hd)) ** 2).sum()
        if best is None or cv < best[0]: best = (cv, h)
    return best


def sigma(z):
    return 1.0 / (1.0 + np.exp(-z))


def softplus(z):
    return np.logaddexp(0.0, z)


def e6_core(x, yA, SA, yB, Ys, truth_d, tag="E6"):
    """x log t window, yA log pA, SA bootstrap cov of log pA, yB, Ys (B x J bootstrap of log pB)."""
    J = len(x)
    SAs = (1 - SHRINK) * SA + SHRINK * np.diag(np.diag(SA)) + 1e-10 * np.trace(SA) / J * np.eye(J)
    SAi = np.linalg.inv(SAs)
    SB = np.cov(Ys, rowvar=False); SBs = SB + 1e-10 * np.trace(SB) / J * np.eye(J)
    SBi = np.linalg.inv(SBs)
    ms = np.linspace(x[0] - 1.0, x[-1] + 1.0, M_GRID_N)
    base = np.stack([np.ones(J), x], 1)
    fits = []
    for m in ms:
        for w in W_GRID:
            phi = 0.5 * w * softplus((x - m) / w)
            X = np.column_stack([base, phi])
            Kx = np.linalg.solve(X.T @ SAi @ X, X.T @ SAi)
            res = yA - X @ (Kx @ yA)
            fits.append((m, w, float(res @ SAi @ res)))
    rss = np.array([f[2] for f in fits]); keep = rss - rss.min() <= CHI2_2_99
    Lo = np.full(J - 1, np.inf); Hi = np.full(J - 1, -np.inf)
    cLo = np.full(J, np.inf); cHi = np.full(J, -np.inf)
    cr = np.sqrt(CHI2_2_96)
    for (m, w, _), k in zip(fits, keep):
        if not k: continue
        phi = 0.5 * w * softplus((x - m) / w)
        X = np.column_stack([base, phi])
        K = np.linalg.solve(X.T @ SBi @ X, X.T @ SBi)[1:]          # rows: a (=-d1/2), H
        th = K @ yB; Sg = K @ SB @ K.T
        s = sigma((x - m) / w)
        # d_j = d1 - H s_j = -2a - H s_j  ; v = (-2, -s_j) on (a,H)
        V = np.stack([np.full(J, -2.0), -s], 1)
        cen = V @ th; sdv = np.sqrt(np.maximum(np.einsum("ji,ik,jk->j", V, Sg, V), 0))
        cLo = np.minimum(cLo, cen - cr * sdv); cHi = np.maximum(cHi, cen + cr * sdv)
        V2 = np.stack([np.zeros(J - 1), -(s[1:] - s[0])], 1)           # Delta_j = -H (s_j - s_0)
        cen2 = V2 @ th; sd2 = np.sqrt(np.maximum(np.einsum("ji,ik,jk->j", V2, Sg, V2), 0))
        Lo = np.minimum(Lo, cen2 - cr * sd2); Hi = np.maximum(Hi, cen2 + cr * sd2)
    r = dict(est=tag, abstain=0, J=J, hw=float(np.mean(cHi - cLo) / 2), n_set=int(keep.sum()), n_grid=len(fits),
             decl_neg=int((Hi < -DELTA).any()), decl_pos=int((Lo > DELTA).any()), decl_zero=int(((Hi < 0) | (Lo > 0)).any()),
             dmin=float(Lo.min()), dmax=float(Hi.max()))
    if truth_d is not None:
        r["cov_simult"] = int(((truth_d >= cLo) & (truth_d <= cHi)).all()); r["cov_point"] = float(((truth_d >= cLo) & (truth_d <= cHi)).mean())
        dt = truth_d[1:] - truth_d[0]
        r["delta_cov"] = int(((dt >= Lo) & (dt <= Hi)).all()); r["true_dmin"] = float(dt.min())
    return r


def est_E5E6(Ublk, mwalk, ts, N, rngB, rngA, truth, thr=L.THR, B=200):
    """returns [E5qu, E5q, E6]."""
    tags = ["E5qu", "E5q", "E6"]
    UA = [u[0::2] for u in Ublk]; UB = [u[1::2] for u in Ublk]
    pA = grand_mean(UA); WA = sum(u.shape[0] for u in UA) * mwalk; cntA = pA * WA
    J = window_len(cntA, pA, N, thr)
    if J < 6:
        return [abst(t, "window") for t in tags]
    x = np.log(ts[:J].astype(float)); yA = np.log(pA[:J]); wA = cntA[:J]
    best = slope_ops(yA, wA, x)
    if best is None:
        return [abst(t, "cvfail") for t in tags]
    hcv = best[1]
    WB = sum(u.shape[0] for u in UB) * mwalk; floor = 0.5 / WB
    pB = np.maximum(grand_mean(UB), floor)[:J]; yB = np.log(pB)
    Ps = np.maximum(boot_means(UB, B, rngB)[:, :J], floor); Ys = np.log(Ps)
    dtrue = truth["d_exact"](ts[:J]) if truth is not None else None
    out = []
    for tag, hh in (("E5qu", 0.7 * hcv), ("E5q", hcv)):
        _, sl = A.poly_rows(x, wA, hh, 2)
        r = e5_core(-2 * sl @ yB, -2 * Ys @ sl.T, dtrue, tag, True); r["h"] = hh; r["tlast"] = int(ts[J - 1]); out.append(r)
    # E6
    PsA = np.maximum(boot_means(UA, B, rngA)[:, :J], 0.5 / WA)
    SA = np.cov(np.log(PsA), rowvar=False)
    r6 = e6_core(x, yA, SA, yB, Ys, dtrue); r6["tlast"] = int(ts[J - 1]); out.append(r6)
    return out
